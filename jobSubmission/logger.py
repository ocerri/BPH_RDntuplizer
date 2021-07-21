import os
import sys
from datetime import datetime

DEBUG = 0
VERBOSE = 1
NOTICE = 2
WARNING = 3

class Logger(object):
    """
    Simple logger class that I wrote for the SNO+ DAQ. Very easy to use:

        log = Logger()
        log.set_logfile("test.log")
        log.notice("blah")
        log.warn("foo")

    The log file format is taken from the Redis log file format which is really
    nice since it shows the exact time and severity of each log message.
    """
    def __init__(self):
        self.logfile = sys.stdout
        # by default, we log everything
        self.verbosity = DEBUG

    def set_verbosity(self, level):
        if isinstance(level, int):
            self.verbosity = level
        elif isinstance(level, str):
            if level == 'debug':
                self.verbosity = DEBUG
            elif level == 'verbose':
                self.verbosity = VERBOSE
            elif level == 'notice':
                self.verbosity = NOTICE
            elif level == 'warning':
                self.verbosity = WARNING
            else:
                raise ValueError("unknown loglevel '%s'" % level)
        else:
            raise TypeError("level must be a string or integer")

    def set_logfile(self, filename):
        self.logfile = open(filename, 'a')

    def debug(self, msg):
        self.log(DEBUG, msg)

    def verbose(self, msg):
        self.log(VERBOSE, msg)

    def notice(self, msg):
        self.log(NOTICE, msg)

    def warn(self, msg):
        self.log(WARNING, msg)

    def log(self, level, msg):
        if level < self.verbosity:
            return

        c = '.-*#'
        pid = os.getpid()
        now = datetime.now()
        buf = now.strftime('%d %b %H:%M:%S.%f')[:-3]

        self.logfile.write('%d:%s %c %s\n' % (pid, buf, c[level], msg))
        self.logfile.flush()

