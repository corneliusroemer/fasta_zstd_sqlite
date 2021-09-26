from setuptools import setup

setup(
    name="fasta_zstd_sqlite",
    description="CLI tool to store and retrieve individual zstd compressed FASTA records in sqlite",
    license="MIT",
    author="Cornelius Roemer",
    url="https://github.com/corneliusroemer/fasta_zstd_sqlite",
    version="0.0.1",
    py_modules=["fasta_zstd_sqlite"],
    install_requires=[
        "Click>=8",
        "Cloup>=0.11",
        "sqlalchemy>=1.3.0",
        "pyzstd>=0.14",
        "biopython>=1.71",
    ],
    python_requires=">=3.7",
    entry_points={
        "console_scripts": [
            "fzs = main:cli",
        ],
    },
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Environment :: Console",
        "Operating System :: POSIX",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.9",
        "Topic :: Utilities",
    ],
    keywords="cli"
)
