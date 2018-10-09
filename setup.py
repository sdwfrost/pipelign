from setuptools import find_packages, setup

with open('README.md', 'r') as f:
    long_description = f.read()

setup(
    name='pipelign',
    version='0.1',
    author='Mukarram Hossain',
    author_email='asmmh2@cam.ac.uk',
    packages=find_packages(exclude=['docs']),
    url='https://github.com/asmmhossain/pipelign',
    license='MIT',
    description='A pipeline for automated alignment',
    long_description=long_description,
    install_requires=[
        'biopython',
        'ete3'
    ],
    scripts=[
        'bin/pipelign'
    ]
)
