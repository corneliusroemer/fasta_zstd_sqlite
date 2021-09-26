import click

import sqlalchemy
from sqlalchemy import create_engine,Table, Column, Integer, String
from sqlalchemy.orm import sessionmaker, declarative_base, relationship, Session
from sqlalchemy.sql.sqltypes import BLOB
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyzstd import ZstdDict, compress, decompress, train_dict, finalize_dict

Base = declarative_base()

class Fasta(Base):
    __tablename__ = 'fasta'

    strain = Column(String, primary_key=True)
    sequence = Column(BLOB)

    def __repr__(self):
        return f"Sequence(strain={self.strain!r}, name={self.sequence!r})"
class Zstd_dict_table(Base):
    __tablename__ = 'zstd_dict'

    id = Column(String, primary_key=True)
    dictionary = Column(BLOB)

    def __repr__(self):
        return f"Dictionary(id={self.id!r}, name={ZstdDict(self.dictionary)!r})"

def vacuum(engine):
    with engine.begin() as conn:
        conn.execute("VACUUM")

def connect_to_db(path, debug=False):
    return create_engine(f"sqlite+pysqlite:///{path}",echo=debug)

def drop_all(engine):
    Base.metadata.drop_all(engine)
    vacuum(engine)

def create_tables(engine):
    Base.metadata.create_all(engine)

def start_session(engine)->Session:
    Session = sessionmaker(bind=engine)
    return Session()

def close_session(session):
    session.close()

def store_dict(session, zd):
    """Store dictionary in table Zstd_dict"""
    zd_dict = Zstd_dict_table(id=zd.dict_id, dictionary=zd.dict_content)
    session.add(zd_dict)
    session.commit()

def add_fasta(session, strain, sequence, zd=None):
    fasta = Fasta(strain=strain, sequence=compress(sequence.encode('UTF-8'),zstd_dict=zd))
    session.add(fasta)

def write_fasta(session, fasta_path, strains=None):
    provided_dict = session.query(Zstd_dict_table).one_or_none()
    if provided_dict:
        zd = ZstdDict(provided_dict.dictionary)
    if strains:
        sequences = session.query(Fasta).filter(Fasta.strain.in_(strains)).all()
    else:
        sequences = session.query(Fasta).all()
    with open(fasta_path, "w") as fp:
        for sequence in sequences:
            record = SeqRecord(
                Seq(decompress(sequence.sequence,zstd_dict=zd).decode('UTF-8')),
                id=sequence.strain,
                description=''
            )
            SeqIO.write(record, fp, "fasta-2line")

@click.group()
def cli():
    pass

@cli.command()
@click.option('--fasta-path', required=True, help='Fasta should be xz compressed')
@click.option('--db-path', required=True)
@click.option('--dict-path')
def store(fasta_path, db_path, dict_path):
    """Program that reads in a fasta and stores it in sqlite, sequence by sequence"""
    #TODO read fasta from stdin
    engine = connect_to_db(db_path)
    drop_all(engine)
    create_tables(engine)
    session = start_session(engine)
    zd = None
    if dict_path:
        with open(dict_path, 'rb') as f:
            file_content = f.read()
        zd = ZstdDict(file_content)
        store_dict(session,zd)
    for record in SeqIO.parse(fasta_path, "fasta"):
        add_fasta(session, record.id, str(record.seq), zd)
    session.commit()
    close_session(session)
    vacuum(engine)

@cli.command()
@click.option('--db-path', required=True)
@click.option('--fasta-path', required=True)
@click.option('--strains-path')
@click.option('--debug', is_flag=True, default=False)
def retrieve(db_path, fasta_path, strains_path, debug):
    """Read in a db and writes out a fasta"""
    engine = connect_to_db(db_path, debug)
    session = start_session(engine)
    # Turn strains file into list of strains
    # Query sequences in list of strains
    # Load strains.txt and turn into list of strains
    strains = None
    if strains_path:
        with open(strains_path, 'r') as f:
            strains = f.read().splitlines()
    write_fasta(session, fasta_path, strains)
    close_session(session)

@cli.command()
@click.option('--db-path', required=True)
@click.option('--dict-path', required=True)
def generate_dict(db_path,dict_path):
    """Create dictionary from provided database"""
    engine = connect_to_db(db_path)
    session = start_session(engine)
    def samples():
        sequences = session.query(Fasta).all()
        for sequence in sequences:
            yield bytes(decompress(sequence.sequence).decode('UTF-8'),'UTF-8')
    dict_size = 100*1024
    raw_dict = train_dict(samples(), dict_size)
    final_dict = finalize_dict(raw_dict, samples(), dict_size, 3)
    with open(dict_path, 'wb') as fp:
        fp.write(final_dict.dict_content)
    close_session(session)

if __name__ == '__main__':
    cli()