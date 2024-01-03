"""
plastic crystal predictor: roughly predict new plastic crystal structures using AIRSS and the machine learned CHGNET.
"""

from setuptools import find_packages, setup


with open("README.md") as file:
    long_description = file.read()

setup(
    name="plastic-crystal-predictor",
    version=0.1,
    description="Plastic-Crystal-Predictor",
    url="https://github.com/badw/plastic-crystal-predictor",
    author="Benjamin A. D. Williamson",
    author_email="benjamin.williamson@ntnu.no",
    long_description=long_description,
    license="MIT",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    keywords="chemistry ase vasp structure",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "chgnet",
        "pyxtal",
        "pymatgen",
        "ase",
        "numpy",
        "smilesbox",
        'pathos',
        'pandas',
        'itertools',
        'airsspy',
        'pickle'
    ],
    #entry_points={
    #    "console_scripts":[
    #        "pcp-run = pcp.cli:main"
    #        ]
    #    }
)
