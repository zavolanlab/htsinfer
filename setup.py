"""HTSinfer package definition"""

from setuptools import setup, find_packages

# Read long description from file
with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

# Read requirements from file
INSTALL_REQUIRES = []
with open("requirements.txt") as fh:
    INSTALL_REQUIRES = fh.read().splitlines()

setup(
    name="htsinfer",
    version="0.1.1",
    description=(
        "HTSinfer infers metadata from High Throughput Sequencing (HTS) data"
    ),
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url="https://github.com/zavolanlab/htsinfer",
    author="Rohan Kandhari",
    author_email="rohan.kandhari.bme16@iitbhu.ac.in",
    maintainer="Rohan Kandhari",
    maintainer_email="rohan.kandhari.bme16@iitbhu.ac.in",
    classifiers=[
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Utilities",
    ],
    entry_points={
        'console_scripts': [
            'htsinfer = src.htsinfer:main',
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
    install_requires=INSTALL_REQUIRES,
    include_package_data=True,
    setup_requires=[
        "setuptools_git == 1.2",
    ],
)
