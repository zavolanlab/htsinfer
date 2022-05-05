"""HTSinfer package definition"""

from setuptools import setup, find_packages

from htsinfer import __version__

# Read long description from file
with open("README.md", "r", encoding="utf-8") as fh:
    LONG_DESCRIPTION = fh.read()

setup(
    name="htsinfer",
    version=__version__,
    description=(
        "Infer experiment metadata from High Throughput Sequencing (HTS) data"
    ),
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url="https://github.com/zavolanlab/htsinfer",
    author="Rohan Kandhari",
    author_email="rohan.kandhari.bme16@iitbhu.ac.in",
    maintainer="Alexander Kanitz",
    maintainer_email="alexander.kanitz@alumni.ethz.ch",
    classifiers=[
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Utilities",
    ],
    entry_points={
        'console_scripts': [
            'htsinfer = htsinfer.cli:main',
        ],
    },
    keywords=[
        'bioinformatics',
        'ngs',
        'high-throughput sequencing',
        'inference',
    ],
    project_urls={
        "Repository": "https://github.com/zavolanlab/htsinfer",
        "Tracker": "https://github.com/zavolanlab/htsinfer/issues",
    },
    packages=find_packages(),
    include_package_data=True,
    setup_requires=[
        "setuptools_git == 1.2",
    ],
)
