#! /usr/bin/env python

import os
import sys
import tempfile
import gzip

from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

def expand_path(path):
    return os.path.abspath(os.path.expandvars(os.path.expanduser(path)))

def is_gzipped(file_path):
    l = ''
    try:
        fs = open(expand_path(file_path), 'r')
        l = fs.next()
    except:
        return False
    finally:
        fs.close()
    if l.startswith("\x1f\x8b"):
        return True
    return False

def process_file_arg(file_arg, mode='rU', compresslevel = None):
    if isinstance(file_arg, str):
        return OpenFile(file_arg, mode, compresslevel), True
    return file_arg, False

class OpenFile(object):
    open_files = set()
    def __init__(self, *args, **kwargs):
        self.file_stream = None
        self.name = None
        self.mode = None
        self.buffering = None
        self.compresslevel = None
        self.gzipped = False
        if kwargs.get('compresslevel', None):
            self.compresslevel = kwargs.pop('compresslevel')
        self.name = kwargs.get('name', None)
        if not self.name:
            self.name = args[0]
        self.name = expand_path(self.name)
        if len(args) > 1:
            self.mode = args[1]
        else:
            self.mode = kwargs.get('mode', 'r')
        if len(args) > 2:
            self.buffering = args[2]
        else:
            self.buffering = kwargs.get('buffering', None)
        if os.path.isfile(self.name):
            self.gzipped = is_gzipped(self.name)

        if self.compresslevel:
            self.mode = self.mode.strip('U+') + 'b'
            self.file_stream = gzip.open(filename = self.name,
                    mode = self.mode,
                    compresslevel = self.compresslevel)
        elif self.gzipped:
            self.mode = self.mode.strip('U+') + 'b'
            self.compresslevel = 9
            self.file_stream = gzip.open(filename = self.name,
                    mode = self.mode,
                    compresslevel = self.compresslevel)
        else:
            if self.buffering:
                self.file_stream = open(name = self.name, mode = self.mode,
                        buffering = self.buffering)
            else:
                self.file_stream = open(name = self.name, mode = self.mode)
        self.__class__.open_files.add(self.name)

    def close(self):
        self.file_stream.close()
        self.__class__.open_files.remove(self.name)

    def __enter__(self):
        return self.file_stream

    def __exit__(self, type, value, traceback):
        self.close()

    def __iter__(self):
        return self.file_stream.__iter__()

    def write(self, string):
        self.file_stream.write(string)

    def writelines(self, sequence_of_strings):
        self.file_stream.writelines(sequence_of_strings)

    def next(self):
        return self.file_stream.next()

    def read(self):
        return self.file_stream.read()

    def readline(self):
        return self.file_stream.readline()

    def readlines(self):
        return self.file_stream.readlines()

    def _get_closed(self):
        return self.file_stream.closed

    closed = property(_get_closed)

    def seek(self, i):
        self.file_stream.seek(i)

    def flush(self):
        self.file_stream.flush()

    def fileno(self):
        return self.file_stream.fileno()

class TemporaryFilePath(object):
    def __enter__(self):
        fd, self.temp_path = tempfile.mkstemp()
        os.close(fd)
        return self.temp_path

    def __exit__(self, type, value, traceback):
        if os.path.exists(self.temp_path):
            os.remove(self.temp_path)

