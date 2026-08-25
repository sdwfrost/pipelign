"""Command-line interface for pipelign."""

'''
    Pipelign.py

    A python based program to align virus sequences.

    The program takes as input a single FASTA formatted file

    Returns a FASTA formatted alignment file

    Different flavours are present:
        - Built-in subprocess concurrency
        - Joblib for Linux, MacOS and Windows
'''

#*********************************************************************

import sys, os, shutil, subprocess, argparse, textwrap, stat, time
from glob import glob
from pipelign import *
from pipelign import _validate_fasta_artifact
from pipelign.reporting import RunReport
from pipelign.skill import skill_main

_active_report = None


def _record_alignment_summary(path, clusters=None, fragments=None):
    records = _validate_fasta_artifact(path, "final alignment", alignment=True)
    details = {
        "final_sequences": len(records),
        "alignment_columns": len(records[0].seq),
    }
    if clusters is not None:
        details["clusters"] = clusters
    if fragments is not None:
        details["fragments"] = fragments
    _active_report.update_summary(**details)


def _artifact_paths(*patterns):
    return sorted(
        {
            path
            for pattern in patterns
            for path in glob(pattern)
            if os.path.isfile(path)
        }
    )


def _pipeline_main(argv=None):
    global _active_report
    configure_resume(False)
    argv = sys.argv[1:] if argv is None else list(argv)
    if argv and argv[0] == "skill":
        return skill_main(argv[1:])

    # get all the command line arguments
    args = getArguments(argv)

    cDir = os.getcwd() # save path for the current working directory

    mArgs = MyStruct(
        inFile = args.inFile,
        outFile = os.path.abspath(args.outFile),
        lenThr = args.lenThr,
        #gCode = args.code,
        alphabet = args.alphabet.lower(),
        keepOrphans = args.keepOrphans,
        keepBadSeqs = args.keepBadSeqs,
        makeZip = args.mZip,
        simPer = args.simPer,
        run = args.run.upper(),
        thread = args.thread,
        mIterL = args.mIterateLong,
        mIterM = args.mIterateMerge,
        fragEmpty = 1,
        longName = 'long.fas',
        fragName = 'frag.fas',
        outDir = args.outDir,
        ambigPer = args.ambigPer,
        merge = args.merge.upper(),
        cFlag = args.clearExistingDirectory,
        stage = args.stage,
        excludeClusters = args.excludeClusters,
        resume = args.resume)

    '''
    # check whether input file exists
    if not checkPresenceOfFile(mArgs.inFile):
        msg = '\n==============================='
        msg += '\nInput file <%s> does not exist in the given path' % mArgs.inFile
        msg += '\nPlease try again with the correct file name and path'
        msg += '\n==============================='
        msg += '\n\nPipelign is exiting.\n'
        sys.exit(msg)
    '''
    valid = validateFastaType(mArgs.inFile,mArgs.alphabet)

    if not valid:
        msg = '\nERROR: Input file does not seem to contain {} residues.'.format(mArgs.alphabet)
        msg += ' Please check you provide correct alphabet type with -a\n'
        sys.exit(msg)

    check_external_dependencies(mArgs.alphabet, mArgs.run, mArgs.stage)

    # Pipeline process starts here
    tFileName = 'input.fas'  # input file in temporary directory
    deFileName = 'input.dealign.fas' # dealigned input file

    '''
    # get current time and generate name for the zipped temporary directory
    timeNow = time.strftime('%Y-%m-%d-%H%M%S')
    zName = 'pipelign.' + timeNow
    '''
    oName = mArgs.outDir

    if mArgs.resume:
        if not os.path.isdir(oName):
            sys.exit("\nERROR: Cannot resume because output directory <%s> does not exist\n" % oName)
        _active_report = RunReport.resume(
            output_dir=os.path.abspath(oName),
            input_path=os.path.abspath(mArgs.inFile),
            output_path=mArgs.outFile,
            arguments=vars(args),
            archive=mArgs.makeZip,
        )
    else:
        makeOutDir(oName,mArgs.cFlag)
        _active_report = RunReport(
            output_dir=os.path.abspath(oName),
            input_path=os.path.abspath(mArgs.inFile),
            output_path=mArgs.outFile,
            arguments=vars(args),
            tools=required_external_tools(mArgs.alphabet, mArgs.run, mArgs.stage),
            archive=mArgs.makeZip,
        )
    configure_resume(
        mArgs.resume,
        _active_report.artifact_checksums() if mArgs.resume else None,
    )

    # create temporary directory
    #tempDir, tDirName = makeTempDir(mArgs.tempDirPath,zName)


    # copy input file inside the output directory
    dName = os.path.join(oName, tFileName)
    copyFile(mArgs.inFile,dName)

    msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Pipelign started with sequence file <%s>\n' % mArgs.inFile
    print(msg)

    with open(os.path.join(oName,'pipelign.log'),'a') as oh:
        oh.write(msg)

    # change current working directory to the output directory
    os.chdir(oName)

    # removes any possible gaps from the sequence file
    deAlign(tFileName, deFileName)

    # make separate sequence file for long sequences and fragments
    mArgs.fragEmpty, numFrag = separateFullFragment(deFileName, mArgs.lenThr, mArgs.longName, mArgs.fragName)

    # separate the long sequences containing ambiguous characters more than a given proportion
    separateLongAmbigs(mArgs.alphabet,mArgs.ambigPer)
    mArgs.longName = 'long.good.fas'

    # use CD-HIT to group sequences based on similarity
    runCDHIT(mArgs.longName, mArgs.alphabet, mArgs.simPer, mArgs.thread)

    # create separate files for each cluster
    numClusters = makeClusters(mArgs.longName,'grp.clstr')
    # add cluster numbers and size to header
    addClusterNumberToReps('grp','long.ClusterList.txt','clsReps.fas')
    _active_report.checkpoint(
        "clustering",
        clusters=numClusters,
        fragments=numFrag,
        artifacts=_artifact_paths(
            "grp",
            "grp.clstr",
            "long.ClusterList.txt",
            "clsReps.fas",
            "long.*.fas",
            "long.good.fas",
        ),
    )


    if numClusters < 3: # less than 3 clusters, tree will not be built
        msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' Only %d cluster(s) created. Phylogenetic tree cannot be built\n' % numClusters
        print(msg)

    else: # at least 3 clusters
        # ********************************************************
        # align cluster reps; reconstruct phylogeny; draw midpoint rooted tree
        msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' Aligning cluster representative sequences\n'
        print(msg)

        if numClusters <= 100:
            linsi('clsReps.fas','clsReps.aln',mArgs.thread,mArgs.mIterL,'pipelign.log')
        else:
            fftnsi('clsReps.fas','clsReps.aln',mArgs.thread,mArgs.mIterL,'pipelign.log')

        msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' Alignment of cluster representatives written in <clsReps.aln>\n'
        print(msg)

        # generate phylogenetic tree using IQ-TREE
        makeIQTree('clsReps.aln',mArgs.thread,mArgs.alphabet)

        msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' Phylogenetic tree of the cluster representatives written in <clsReps.aln.treefile>\n'
        print(msg)

        # make midpoint root tree for cluster representatives
        makeMidPointRootTree('clsReps.aln.treefile')

        msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' Midpoint rooted tree written in <clsReps.aln.midpoint.treefile>\n'
        print(msg)

        # draw midpoint root tree
        print('\nMidpoint rooted tree for cluster representatives is drawn below. Tips are labelled as <H_N_C>.')
        print('\tH: header of the cluster representative sequence')
        print('\tN: number of long sequences in that cluster')
        print('\tC: unique identifier of the cluster')
        drawAsciiTree('clsReps.aln.midpoint.treefile')
        print('\n')
        # ********************************************************

    # make individual long cluster alignments
    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Aligning long sequence clusters\n'
    print(msg)

    # align long sequence clusters and create <long.xx.aln> files
    if mArgs.run == 'G': # run built-in concurrent version
        alnLongSequenceClustersGNUParallel(numClusters,mArgs.thread,mArgs.mIterL)

    elif mArgs.run == 'J':
        alnLongSequenceClustersJoblibParallel(numClusters,mArgs.thread,mArgs.mIterL)

    # make HMMs and database from long sequence cluster alignments
    makeHMMdbParallel(numClusters,mArgs.thread,mArgs.alphabet,'hmmer.log')
    _active_report.checkpoint(
        "cluster_alignments_and_hmms",
        clusters=numClusters,
        artifacts=_artifact_paths(
            "long.*.aln",
            "long.*.hmm",
            "pipelign.hmm*",
            "clsReps.aln*",
        ),
    )

    msg = '\nAlignment files and HMMs can be found in <%s>\n' % oName
    msg += '\tLong sequence alignments have names <long.xx.aln>\n'
    msg += '\tHMM file written in <long.xx.hmm>\n'
    msg += '\tHMM database written in <pipelign.hmm>\n'

    print(msg)
    with open('pipelign.log','a') as oh:
        oh.write(msg)
    #***************
    # only create long sequence cluster alignments and HMMs and exit
    if mArgs.stage == 1:
        print('\nThank you for using Pipelign.\n')
        sys.exit(0)

    #***************
    # create alignment from only long sequences
    if mArgs.stage == 2:
        if mArgs.merge == 'P':
            longSeqAlignmentParallel(numClusters,mArgs)

        elif mArgs.merge == 'C':
            longSeqAlignmentConsensus(numClusters,mArgs)

        # remove clusters
        if mArgs.excludeClusters:
            excludeClustersFromalignment('final.all.aln','final.aln')
        else:
            copyFile('final.all.aln','final.aln')

        if not checkPresenceOfFile('final.aln'):
            raise PipelignError('Stage 2 did not produce <final.aln>')
        publishAlignment('final.aln',mArgs.outFile)
        _record_alignment_summary('final.aln', numClusters, numFrag)

        msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' Final alignment written in <%s>\n' % mArgs.outFile
        print(msg)

        print('\nThank you for using Pipelign.\n')
        sys.exit(0)

    #***************
    # Assigns clusters to fragments based on BLAST and HMM search
    if mArgs.stage == 3:
        msg = ''
        # if no fragments present in the dataset
        if not checkPresenceOfFile('frag.fas'):
            msg += '\nFragment file <frag.fas> does not exist. No fragments found for assigning clusters.\n'
            print(msg)
            with open('pipelign.log', 'a') as oh:
                oh.write(msg)

        else:
            assignClusterstoFrags(numClusters,numFrag,mArgs)
            _active_report.checkpoint(
                "fragment_assignment",
                fragments=numFrag,
                artifacts=_artifact_paths(
                    "frag.ClusterList.txt",
                    "frag.*.fas",
                    "frag.blast.txt",
                    "hmm.out",
                ),
            )

            msg += '\tCluster files for fragments are written in <frag.xx.fas>\n'

        print('\nThank you for using Pipelign.\n')
        sys.exit(0)


    #***************
    # make cluster alignments with fragments
    if mArgs.stage == 4:
        msg = ''
        # if no fragments present in the dataset
        if not checkPresenceOfFile('frag.fas'):
            msg += '\nFragment file <frag.fas> does not exist. No fragments found for assigning clusters.'
            print(msg)
            with open('pipelign.log', 'a') as oh:
                oh.write(msg)


        else:
            assignClusterstoFrags(numClusters,numFrag,mArgs)
            _active_report.checkpoint(
                "fragment_assignment",
                fragments=numFrag,
                artifacts=_artifact_paths(
                    "frag.ClusterList.txt",
                    "frag.*.fas",
                    "frag.blast.txt",
                    "hmm.out",
                ),
            )

            if mArgs.run == 'G':
                addFragmentToClustersGNUParallel(numClusters,'fragAlign.log',mArgs)
            elif mArgs.run == 'J':
                addFragmentToClustersJoblibParallel(numClusters,'fragAlign.log',mArgs.thread,mArgs.mIterL)

            msg = '\tCluster files for fragments are written in <frag.xx.fas>\n'

            msg += '\tCluster alignments are written in <cls.xx.aln>\n'

            print(msg)
            with open('pipelign.log', 'a') as oh:
                oh.write(msg)


        print('\nThank you for using Pipelign.\n')
        sys.exit(0)


    #************
    # complete alignment requested, -n 5
    if mArgs.stage == 5:
        msg = ''
        # if no fragments present in the dataset
        if not checkPresenceOfFile('frag.fas'):
            msg += '\nFragment file <frag.fas> does not exist. No fragments found for assigning clusters.'
            print(msg)
            with open('pipelign.log', 'a') as oh:
                oh.write(msg)

        else:
            assignClusterstoFrags(numClusters,numFrag,mArgs)
            _active_report.checkpoint(
                "fragment_assignment",
                fragments=numFrag,
                artifacts=_artifact_paths(
                    "frag.ClusterList.txt",
                    "frag.*.fas",
                    "frag.blast.txt",
                    "hmm.out",
                ),
            )
            msg = 'Cluster files for fragments are written in <frag.xx.fas>\n'
            print(msg)
            with open('pipelign.log', 'a') as oh:
                oh.write(msg)


        if mArgs.run == 'G':
            addFragmentToClustersGNUParallel(numClusters,'fragAlign.log',mArgs)
        elif mArgs.run == 'J':
            addFragmentToClustersJoblibParallel(numClusters,'fragAlign.log',mArgs.thread,mArgs.mIterL)
        _active_report.checkpoint(
            "fragment_alignments",
            clusters=numClusters,
            artifacts=_artifact_paths("cls.*.aln"),
        )

        msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' Cluster alignments are written in <cls.xx.aln>\n'
        print(msg)
        with open('pipelign.log', 'a') as oh:
            oh.write(msg)


        if numClusters == 1:
            copyFile('cls.0.aln','final.noOrphans.aln')
        #elif numClusters > 1:
            #mergeConsensus(numClusters,'final.noOrphans.aln','all',mArgs,'consensus.log',cDir,tFileName,zName)

        # only two clusters, merge the pair
        elif numClusters == 2:
            if mArgs.merge == 'P':
                mergePair('cls.0.aln','cls.1.aln','final.noOrphans.aln',mArgs.thread,mArgs.mIterM,'merge.log')
            elif mArgs.merge == 'C':
                mergeConsensus(numClusters,'final.noOrphans.aln','all',mArgs,'consensus.log')

        # more than 2 clusters
        elif numClusters >= 3:
            if mArgs.merge == 'P':
                resName = merge_gnu_parallel('clsReps.aln.midpoint.treefile','all',mArgs.thread,mArgs.mIterL,'merge.log')
                copyFile(resName,'final.noOrphans.aln')
            elif mArgs.merge == 'C':
                mergeConsensus(numClusters,'final.noOrphans.aln','all',mArgs,'consensus.log')

        msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' Alignment with long sequences and fragments written in <final.noOrphans.aln>\n'
        print(msg)
        with open('pipelign.log','a') as oh:
            oh.write(msg)

        #'''

        # add unclustered fragments
        if mArgs.keepOrphans and checkPresenceOfFile('frag.noClusters.fas'):
            # create alignment by merging cluster alignments
            msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
            msg += ' Adding orphan fragments to final alignment\n'
            print(msg)
            with open('pipelign.log','a') as oh:
                oh.write(msg)

            addFragments('frag.noClusters.fas','final.noOrphans.aln','final.noBadSeq.aln',mArgs.thread,mArgs.mIterL,'merge.log')

            msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
            msg += ' Alignment with unclustered fragments written in <final.noBadSeq.aln>\n'
            print(msg)
            with open('pipelign.log','a') as oh:
                oh.write(msg)

        else:
            copyFile('final.noOrphans.aln','final.noBadSeq.aln')

            msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
            msg += ' <final.noBadSeq.aln> file created\n'
            print(msg)
            with open('pipelign.log','a') as oh:
                oh.write(msg)

        # check if <long.bad.fas> needs to be added to the final alignment
        if mArgs.keepBadSeqs and checkPresenceOfFile('long.bad.fas'):
            addLongSeqs('long.bad.fas','final.noBadSeq.aln','final.all.aln',mArgs.thread,mArgs.mIterL,'addBadSeq.log')
        else:
            copyFile('final.noBadSeq.aln','final.all.aln')

        # if user wants to exclude clusters from the final alignment
        # first show the cluster representative tree
        # then asks for list of clusters separated by space
        # remove all the sequences from the <final.all.aln> that belong to the clusters
        if mArgs.excludeClusters:
            excludeClustersFromalignment('final.all.aln','final.aln')
        else:
            copyFile('final.all.aln','final.aln')

        if checkPresenceOfFile('final.aln'):
            publishAlignment('final.aln',mArgs.outFile)

        if checkPresenceOfFile(mArgs.outFile):
            msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
            msg += ' Final alignment written in <%s>\n' % mArgs.outFile
            print(msg)

            _active_report.checkpoint(
                "final_alignment",
                output=mArgs.outFile,
                artifacts=_artifact_paths(
                    "final.aln",
                    "final.all.aln",
                    "final.noOrphans.aln",
                    "final.noBadSeq.aln",
                ),
            )
            _record_alignment_summary('final.aln', numClusters, numFrag)

            # change to parent directory
            os.chdir(cDir)

        else:
            msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
            msg += ' Pipelign could not write the file <%s>\n' % mArgs.outFile
            msg += '\tSomething went wrong. Please check the output directory.\n'
            print(msg)
            sys.exit(1)
            #cZip(cDir,tFileName,zName)


        print('\nThank you for using Pipelign.\n')
        sys.exit(0)


def main(argv=None):
    global _active_report
    _active_report = None
    original_cwd = os.getcwd()
    try:
        return _pipeline_main(argv)
    except SystemExit as exc:
        if _active_report is not None:
            code = exc.code if isinstance(exc.code, int) else 1
            try:
                _active_report.finish(
                    "success" if code == 0 else "failed",
                    None if code == 0 else str(exc),
                )
            except Exception as finish_exc:
                print("pipelign: error: %s" % finish_exc, file=sys.stderr)
                raise SystemExit(1) from finish_exc
        raise
    except Exception as exc:
        if _active_report is not None:
            _active_report.finish("failed", "%s: %s" % (type(exc).__name__, exc))
        print("pipelign: error: %s" % exc, file=sys.stderr)
        return 1
    finally:
        os.chdir(original_cwd)


if __name__=="__main__":
    raise SystemExit(main())
