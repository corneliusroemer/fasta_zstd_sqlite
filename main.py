from io import TextIOWrapper
import click
from sqlalchemy import create_engine,Table, Column, Integer, String
from sqlalchemy.orm import sessionmaker, declarative_base, relationship, Session
from sqlalchemy.sql.sqltypes import BLOB
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyzstd import ZstdDict, compress, compressionLevel_values, decompress, train_dict, finalize_dict
import cloup
from cloup import group, command, option, option_group
from cloup.constraints import require_one
import sys
import contextlib

COMP_LEVEL = 0
DICT_SIZE = 128*1024
DICT_SAMPLE_NUMBER = 1000

@contextlib.contextmanager
def smart_open(filename=None, mode='w'):
    if filename and filename != '-':
        fh = open(filename, 'w')
    else:
        fh = sys.stdout

    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()

@contextlib.contextmanager
def smart_in(filename=None, mode='r'):
    if type(filename) != TextIOWrapper:
        fh = open(filename, 'r')
    else:
        fh = filename

    try:
        yield fh
    finally:
        if type(filename) != TextIOWrapper:
            fh.close()

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

def load_dict(session,dict_path) -> ZstdDict:
    with open(dict_path, 'rb') as f:
        file_content = f.read()
    zd = ZstdDict(file_content)
    store_dict(session,zd)
    return zd

def add_fasta(session, fasta_path, dict_path=None, zd=None, level=COMP_LEVEL, sample_number=None):  
    if not fasta_path or fasta_path == '-':
        fasta_path = click.get_text_stream('stdin')
    if dict_path is not None:
        zd = load_dict(dict_path)
    for count, record in enumerate(SeqIO.parse(fasta_path, "fasta")):
        if sample_number is not None and count == sample_number:
            gen_dict(level=level,size=DICT_SIZE,number=sample_number,passed_session=session)
            provided_dict = session.query(Zstd_dict_table).one_or_none()
            zd = ZstdDict(provided_dict.dictionary)
            for sequence in session.query(Fasta):
                sequence.sequence = compress(decompress(sequence.sequence),level,zstd_dict=zd)
            session.commit()
            
        if zd is None:
            fasta = Fasta(strain=record.id, sequence=compress(str(record.seq).encode('UTF-8'),level))
        else:
            fasta = Fasta(strain=record.id, sequence=compress(str(record.seq).encode('UTF-8'),level,zstd_dict=zd))
        
        if count % 1000 == 0:
            session.commit()
            print(f"Sequence count: {count}")
        
        session.add(fasta)

def write_fasta(session, fasta_path, strains=None):
    provided_dict = session.query(Zstd_dict_table).one_or_none()
    zd = None
    if provided_dict:
        zd = ZstdDict(provided_dict.dictionary)
    if strains:
        sequences = session.query(Fasta).filter(Fasta.strain.in_(strains)).all()
    else:
        sequences = session.query(Fasta).all()
    with smart_open(fasta_path, "w") as fp:
        for sequence in sequences:
            record = SeqRecord(
                Seq(decompress(sequence.sequence,zstd_dict=zd).decode('UTF-8')),
                id=sequence.strain,
                description=''
            )
            SeqIO.write(record, fp, "fasta-2line")

@cloup.group()
def cli():
    pass

@cli.command()
@cloup.option('--fasta-path')
@cloup.option('--db-path', required=True)
@cloup.option('--dict-path')
@cloup.option('--level',default=COMP_LEVEL)
def insert(fasta_path, db_path, dict_path,level):
    """Program that reads in a fasta and stores it in sqlite, sequence by sequence"""
    #TODO give option to autogenerate dictionary
    engine = connect_to_db(db_path)
    drop_all(engine)
    create_tables(engine)
    session = start_session(engine)
    add_fasta(session=session, dict_path=dict_path, fasta_path=fasta_path, level=level, sample_number=1000)
    session.commit()
    close_session(session)
    vacuum(engine)

@cli.command()
@cloup.option('--db-path', required=True)
@cloup.option('--fasta-path',help='Optional, if not provided, will print to stdout')
@cloup.option('--strains-path',help='Optional, if not provided, will read from stdin')
@cloup.option('--debug', is_flag=True, default=False)
def query(db_path, fasta_path, strains_path, debug):
    """Read in a db and writes out a fasta"""
    engine = connect_to_db(db_path, debug)
    session = start_session(engine)
    strains = None
    if strains_path == '-':
        strains_path = click.get_text_stream('stdin')
    if strains_path:
        with smart_in(strains_path, 'r') as f:
            strains = f.read().splitlines()
    write_fasta(session, fasta_path, strains)
    close_session(session)

def gen_dict(db_path=None,dict_path=None,level=COMP_LEVEL,size=DICT_SIZE,number=DICT_SAMPLE_NUMBER,passed_session=None):    
    """Create dictionary from provided database\nRequires uncompressed sequences"""
    if not passed_session:
        engine = connect_to_db(db_path)
        session = start_session(engine)
    else:
        session = passed_session

    provided_dict = session.query(Zstd_dict_table).one_or_none()
    zd = None
    if provided_dict:
        zd = ZstdDict(provided_dict.dictionary)
    def samples():
        count = 0
        sequences = session.query(Fasta).all()
        for sequence in sequences:
            if count < number:
                yield bytes(decompress(sequence.sequence,zstd_dict=zd).decode('UTF-8'),'UTF-8')
            count += 1
    raw_dict = train_dict(samples(), size)
    final_dict = finalize_dict(raw_dict, samples(), size, level)
    print(final_dict)

    # if no dict in db: store it
    if not provided_dict:
        store_dict(session, final_dict)
        session.commit()

    if not passed_session:
        with open(dict_path, 'wb') as fp:
            fp.write(final_dict.dict_content)
        close_session(session)

@cli.command()
@cloup.option('--db-path', required=True)
@cloup.option('--dict-path', required=True)
@cloup.option('--level', default=COMP_LEVEL)
@cloup.option('--size', default=DICT_SIZE)
@cloup.option('--number', default=DICT_SAMPLE_NUMBER)
def generate_dict(db_path=None,dict_path=None,level=COMP_LEVEL,size=DICT_SIZE,number=DICT_SAMPLE_NUMBER,passed_session=None):
    gen_dict(db_path=db_path,dict_path=dict_path,level=level,size=size,number=number,passed_session=passed_session)


if __name__ == '__main__':
    cli()