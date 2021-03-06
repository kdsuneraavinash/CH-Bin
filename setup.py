# -*- coding: utf-8 -*-
from setuptools import setup

packages = ["ch_bin", "ch_bin.cli", "ch_bin.core", "ch_bin.core.clustering", "ch_bin.core.features"]

package_data = {"": ["*"]}

install_requires = [
    "numpy>=1.21.1,<2.0.0",
    "scipy>=1.7.1,<2.0.0",
    "sklearn>=0.0,<0.1",
    "biopython>=1.79,<2.0",
    "pandas>=1.3.1,<2.0.0",
    "cvxopt>=1.2.6,<2.0.0",
    "quadprog>=0.1.8,<0.2.0",
    "numba>=0.53.1,<0.54.0",
    "tqdm>=4.62.0,<5.0.0",
    "click>=8.0.1,<9.0.0",
]

with open("README.md") as readme_file:
    readme = readme_file.read()

setup(
    author="K. D. Sunera Avinash Chandrasiri",
    author_email="kdsuneraavinash@gmail.com",
    python_requires=">=3.8,<3.10",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    description="Taxonomy Independent Hybrid Metagenomic Binning tool utilizing Convex Hull distance metric",
    entry_points={
        "console_scripts": [
            "ch_bin=ch_bin.ch_bin:run",
        ],
    },
    install_requires=install_requires,
    license="MIT license",
    long_description=readme,
    packages=packages,
    package_data=package_data,
    keywords="ch_bin",
    name="ch_bin",
    url="https://github.com/kdsuneraavinash/CH-Bin",
    version="0.0.3",
    zip_safe=False,
)
