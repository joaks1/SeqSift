#! /usr/bin/env python

import os
import sys
import subprocess
import traceback
from cStringIO import StringIO

from seqsift.utils.errors import WorkerExecutionError
from seqsift.utils import functions
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)


class Worker(object):
    total = 0
    def __init__(self,
            stdout_path = None,
            stderr_path = None,
            append_stdout = False,
            subprocess_kwargs = {},
            expected_exit_code = 0,
            tag = None):
        self.__class__.total += 1
        self.stdout_path = stdout_path
        self.stderr_path = stderr_path
        self.cmd = []
        self.process = None
        self.finished = False
        self.exit_code = None
        self.stdout = None
        self.stderr = None
        self.append_stdout = append_stdout
        self.error = None
        self.trace_back = None
        self.subprocess_kwargs = subprocess_kwargs
        self.expected_exit_code = expected_exit_code
        self.tag = tag

    def get_stderr(self):
        if not self.stderr_path:
            return self.stderr
        try:
            se = open(self.stderr_path, 'rU')
        except IOError, e:
            _LOG.error('Could not open stderr file')
            raise e
        msg = se.read()
        se.close()
        return msg

    def start(self):
        try:
            self._start()
        except Exception as e:
            self.error = e
            f = StringIO()
            traceback.print_exc(file=f)
            self.trace_back = f.getvalue()

    def _start(self):
        try:
            self._pre_process()
        except:
            e = StringIO()
            traceback.print_exc(file=e)
            _LOG.error('Error during pre-processing:\n{0}'.format(
                    e.getvalue()))
            raise
        _LOG.debug('Starting process with following command:\n\t'
                '{0}'.format(' '.join([str(x) for x in self.cmd])))
        if self.stdout_path:
            if self.append_stdout:
                sout = open(self.stdout_path, 'a')
            else:
                sout = open(self.stdout_path, 'w')
        else:
            sout = subprocess.PIPE
        if self.stderr_path:
            serr = open(self.stderr_path, 'w')
        else:
            serr = subprocess.PIPE
        self.process = subprocess.Popen(self.cmd,
                stdout = sout,
                stderr = serr,
                shell = False,
                **self.subprocess_kwargs)
        self.stdout, self.stderr = self.process.communicate()
        self.exit_code = self.process.wait()
        if hasattr(sout, 'close'):
            sout.close()
        if hasattr(serr, 'close'):
            serr.close()
        if ((self.expected_exit_code != None) and
                (self.exit_code != self.expected_exit_code)):
            _LOG.error('execution failed')
            raise errors.WorkerExecutionError('{0} failed.\ninvocation:\n{1}\n'
                    'stderr:\n{2}'.format(
                    self.name, ' '.join([str(x) for x in self.cmd]),
                    self.get_stderr()))
        try:
            self._post_process()
        except:
            e = StringIO()
            traceback.print_exc(file=e)
            _LOG.error('Error during post-processing:\n{0}'.format(
                    e.getvalue()))
            raise
        self.finished = True

    def _post_process(self):
        pass

    def _pre_process(self):
        pass

class GenericWorker(Worker):
    count = 0
    def __init__(self,
            exe,
            args,
            stdout_path = None,
            stderr_path = None,
            append_stdout = False,
            subprocess_kwargs = {},
            expected_exit_code = 0,
            tag = None):
        Worker.__init__(self,
                stdout_path = stdout_path,
                stderr_path = stderr_path,
                append_stdout = append_stdout,
                subprocess_kwargs = subprocess_kwargs,
                expected_exit_code = expected_exit_code,
                tag = tag)
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.exe = functions.get_external_tool(exe)
        if isinstance(args, str):
            args = args.strip().split()
        self.args = [str(x) for x in args]
        self.cmd = [self.exe] + self.args

class PaupWorker(GenericWorker):
    count = 0
    def __init__(self,
            nex_path,
            exe = None,
            stdout_path = None,
            stderr_path = None,
            expected_exit_code = None,
            tag = None):
        self.nex_path = nex_path
        if not exe:
            exe = 'paup'
        GenericWorker.__init__(self,
                exe = exe,
                args = ['-n', self.nex_path],
                stdout_path = stdout_path,
                stderr_path = stderr_path,
                append_stdout = append_stdout,
                subprocess_kwargs = {'shell': False},
                expected_exit_code = expected_exit_code,
                tag = tag)
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)

