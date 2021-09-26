# 1. Specify the UI/API
# 2. Implement API step by step
# 3. Use SQLAlchemy to connect to the database

import click

from sqlalchemy import create_engine
from sqlalchemy import MetaData
from sqlalchemy.orm import sessionmaker
from sqlalchemy import Table, Column, Integer, String
from sqlalchemy.orm import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy.orm import Session
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import sqlalchemy

from pyzstd import compress, decompress

Base = declarative_base()

class Fasta(Base):
    __tablename__ = 'fasta'

    strain = Column(String, primary_key=True)
    sequence = Column(String)

    def __repr__(self):
        return f"Sequence(strain={self.strain!r}, name={self.sequence!r})"

def connect_to_db(path):
    return create_engine(f"sqlite+pysqlite:///{path}", echo=True)

def drop_all(engine):
    Base.metadata.drop_all(engine, checkfirst=True)

def create_tables(engine):
    Base.metadata.create_all(engine)

def start_session(engine)->Session:
    Session = sessionmaker(bind=engine)
    return Session()

def add_fasta(session, strain, sequence):
    #TODO encode sequence with zstd
    fasta = Fasta(strain=strain, sequence=compress(sequence.encode('ascii')))
    session.add(fasta)
    session.commit()

def close_session(session):
    session.close()

def write_fasta(session, fasta_path):
    sequences = session.query(Fasta).all()
    with open(fasta_path, "w") as fp:
        for sequence in sequences:
            record = SeqRecord(
                Seq(decompress(sequence.sequence).decode('ascii')),
                id=sequence.strain,
                description=''
            )
            SeqIO.write(record, fp, "fasta-2line")
            print(sequence.strain)

@click.group()
def cli():
    pass

@cli.command()
@click.option('--fasta-path', required=True, help='Fasta should be xz compressed')
@click.option('--db-path', required=True)
def store(fasta_path, db_path):
    """Program that reads in a fasta and stores it in sqlite, sequence by sequence"""
    engine = connect_to_db(db_path)
    drop_all(engine)
    create_tables(engine)
    session = start_session(engine)
    for record in SeqIO.parse(fasta_path, "fasta"):
        add_fasta(session, record.id, str(record.seq))
    close_session(session)

@cli.command()
@click.option('--db-path', required=True)
@click.option('--fasta-path', required=True)
def retrieve(db_path, fasta_path):
    """Read in a db and writes out a fasta"""
    engine = connect_to_db(db_path)
    session = start_session(engine)
    write_fasta(session, fasta_path)
    close_session(session)

if __name__ == '__main__':
    cli()