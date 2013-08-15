#! /usr/bin/env python

import os
import sys
from gzip import GzipFile

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

class OpenFile(object):
    open_files = set()
    def __init__(self, *args, **kwargs):
        self.file_stream = None
        self.name = None
        self.mode = None
        self.buffering = None
        self.compresslevel = None
        self.gzipped = False
        self.close = True
        if kwargs.get('compresslevel', None):
            self.compresslevel = kwargs.pop('compresslevel')
        self.name = kwargs.get('name', None)
        if not self.name:
            self.name = args[0]
        if isinstance(self.name, str):
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
        else:
            self.file_stream = self.name
            self.close = False

    def __enter__(self):
        return self.file_stream

    def __exit__(self, type, value, traceback):
        if self.close:
            self.file_stream.close()
            self.__class__.open_files.remove(self.name)

