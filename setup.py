#!/usr/bin/env python3
"""
Setup script for Organoverse - AI-assisted workbench for plant organellar genome assembly
"""

from setuptools import setup, find_packages
import os

# Read the contents of README file
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Read requirements
with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name="organoverse",
    version="1.0.0",
    author="Organoverse Development Team",
    author_email="contact@organoverse.org",
    description="AI-assisted workbench for plant mitochondrial and chloroplast genome assembly",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/organoverse/organoverse-workbench",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "organoverse=organoverse.cli:main",
        ],
    },
    include_package_data=True,
    package_data={
        "organoverse": ["data/*", "configs/*"],
    },
)