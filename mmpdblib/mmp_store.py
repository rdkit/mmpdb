# std imports
import json
import unittest
import tempfile
import sys
import os

# local
import sqlitedict
from sqlitedict import SqliteDict

fname = '/home/oriol/dev/mmpdb/tests_dict/d1.sqlite'
db = SqliteDict(filename=fname)

def insert():
    with db:
        db['key'] = 'value'
        db.commit()
    with db:
        db['key2'] = 'value2'
        db.commit()
        db.close()

def read():
    mydict = SqliteDict(fname, autocommit=True)
    return mydict



