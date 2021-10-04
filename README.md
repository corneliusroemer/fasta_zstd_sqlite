CLI tool to store and retrieve individual dictionary zstd compressed FASTA records in an sqlite database

## Features

- Fast retrieval of FASTA records from database (~GB/s)
- High compression ratio (100 GB -> 1 GB) for similar sequences (like Sars-CoV-2 genomes)
- Fast compression (~GB/s)
- Convenient CLI user interface (no need to use Python)
- User can choose between explicit file names or stdin/stdout
- Batteries included:
  - zstd dictionary is auto-generated from first fasta record
  - zstd dictionary is auto-included in sqlite database, read automagically
- sqlite database can be read by external tools
- Records can be queried passing strain names via stdin or strains.txt file

## Roadmap

- Add support to store non-Fasta lines, e.g. metadata rows
- Extend querying capabilities beyond strain names:
  - by metadata
  - by mutations
- Support random sampling
- Store hashes of uncompressed records to allow fast checking which records have changed (save time for downstreaming processing if records are unchanged, e.g. in Sars-CoV-2 pipelines
- Make conda installable
- Add tests

## Installation

1. Clone the repository:
```bash
git clone https://github.com/corneliusroemer/fasta_zstd_sqlite.git
cd fasta_zstd_sqlite
```
2. Install using pip
```python
python3 -m venv env
source env/bin/activate
python3 -m pip install --editable .
```

## Usage

Compressing fasta records line by line in sqlite db:
```
xzcat in.fasta.xz | fzs insert --db-path in.fasta.db
fzs insert --fasta-path in.fasta --db-path in.fasta.db
head in.fasta | fzs insert --db-path in.fasta.db
```

Querying records from sqlite db:
```
fzs query --db-path in.fasta.db --strains-path strains.txt  --fasta-path out.fasta
echo Wuhan-Hu-1 | fzs query --strains-path - --db-path in.fasta.db | less
```

Uncompressing all records back to fasta:
```
fzs query --db-path in.fasta.db --fasta-path out.fasta
```

