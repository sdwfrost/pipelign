# Change Log

All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).
This change log follows the conventions of
[keep a CHANGELOG](http://keepachangelog.com/).

## [Unreleased]

### Added

- Atomic JSON run reports with tool metadata, artifact checksums, summaries,
  checkpoint-compatible resume, and optional ZIP archives.
- Fake-tool integration coverage for stage 1 and stage 5 in both executors.

### Changed

- Replaced generated GNU Parallel shell scripts with bounded Python-managed
  subprocess concurrency; `-r G` remains the compatibility name for this mode.
- Final alignments and archives are validated and atomically published.
- Joblib uses threads for external subprocess tasks so workers share the active
  workspace and resume state.

### Fixed

- RNA and single-record FASTA handling, duplicate identifier detection, HMM
  fallback assignments, HMM single-record aliasing, and malformed tool output.
- External command failures, stale/corrupt resume artifacts, excessive worker
  counts, unsafe permissions, and misleading success reports.

[Unreleased]: https://github.com/asmmhossain/pipelign/compare/v0.1...HEAD
