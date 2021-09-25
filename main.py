# 1. Specify the UI/API
# 2. Implement API step by step
# 3. Use SQLAlchemy to connect to the database

#%%
import click

from sqlalchemy import create_engine
from sqlalchemy import MetaData
from sqlalchemy import Table, Column, Integer, String
from sqlalchemy.orm import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy.orm import Session
from Bio import SeqIO

from pyzstd import compress, decompress

#%%
Base = declarative_base()

class Fasta(Base):
    __tablename__ = 'fasta'

    strain = Column(String, primary_key=True)
    sequence = Column(String)

    def __repr__(self):
        return f"Sequence(strain={self.strain!r}, name={self.sequence!r})"

#%%
def connect_to_db(path) -> Session:
    engine = create_engine(f"sqlite+pysqlite:///{path}", echo=True)
    Base.metadata.drop_all(engine, checkfirst=True)
    Base.metadata.create_all(engine)
    return Session(engine)

def add_fasta(session, strain, sequence):
    #TODO encode sequence with zstd
    fasta = Fasta(strain=strain, sequence=compress(sequence.encode('ascii')))
    session.add(fasta)
    session.commit()

def close_session(session):
    session.close()
#%%
@click.command()
@click.option('--fasta-path', required=True, help='Fasta should be xz compressed')
@click.option('--out-path', required=True)
def store(fasta_path, out_path):
    """Program that reads in a fasta and stores it in sqlite, sequence by sequence"""
    session = connect_to_db(out_path)
    for record in SeqIO.parse(fasta_path, "fasta"):
        add_fasta(session, record.id, str(record.seq))
    close_session(session)

if __name__ == '__main__':
    store()