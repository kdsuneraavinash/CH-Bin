# -*- coding: utf-8 -*-
from setuptools import setup

packages = ["bin_x", "bin_x.analysis", "bin_x.cli", "bin_x.core", "bin_x.core.clustering", "bin_x.core.features"]

package_data = {"": ["*"]}

install_requires = [
    "biopython>=1.79,<2.0",
    "click>=8.0.1,<9.0.0",
    "cvxopt>=1.2.6,<2.0.0",
    "matplotlib>=3.4.2,<4.0.0",
    "numba>=0.53.1,<0.54.0",
    "numpy>=1.21.1,<2.0.0",
    "pandas>=1.3.1,<2.0.0",
    "quadprog>=0.1.8,<0.2.0",
    "scipy>=1.7.1,<2.0.0",
    "seaborn>=0.11.1,<0.12.0",
    "sklearn>=0.0,<0.1",
    "tqdm>=4.62.0,<5.0.0",
]

with open("README.md") as readme_file:
    readme = readme_file.read()

setup(
    author="K. D. Sunera Avinash Chandrasiri",
    author_email="kdsuneraavinash@gmail.com",
    python_requires=">=3.9,<3.10",
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
            "bin_x=bin_x.bin_x:cli",
        ],
    },
    install_requires=install_requires,
    license="MIT license",
    long_description=readme,
    packages=packages,
    package_data=package_data,
    keywords="bin_x",
    name="bin_x",
    url="https://github.com/kdsuneraavinash/bin_x",
    version="0.0.2",
    zip_safe=False,
)
