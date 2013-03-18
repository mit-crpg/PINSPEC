import numpy
from pinspec import *

# logging module that connects with Will's C++
# log.h/log.cpp module
def py_printf(level, my_str, *args):

    level_val = assignValue(level)
    
    if level_val <= get_loglevel:
        if level == 'DEBUG':
            print '[  DEBUG  ]  ' + my_str % args
        elif level == 'INFO':
            print '[  INFO   ]  ' + my_str % args
        elif level == 'NORMAL':
            print '[  NORMAL ]  ' + my_str % args
        elif level == 'WARNING':
            print '[ WARNING ]  ' + my_str % args
        elif level == 'CRITICAL':
            print '[ CRITICAL]  ' + my_str % args
        elif level == 'RESULT':
            print '[  RESULT ]  ' + my_str % args
        elif level == 'ERROR':
            print '[  ERROR  ]  ' + my_str % args



def assignValue(level):

    if level == 'DEBUG':
        return 0
    elif level == 'INFO':
        return 1
    elif level == 'NORMAL':
        return 2
    elif level == 'WARNING':
        return 3
    elif level == 'CRITICAL':
        return 4
    elif level == 'RESULT':
        return 5
    elif level == 'ERROR':
        return 6


def setLevel(level):
    
    if level == 'DEBUG':
        log_setlevel(DEBUG)
    elif level == 'INFO':
        log_setlevel(INFO)
    elif level == 'NORMAL':
        log_setlevel(NORMAL)
    elif level == 'WARNING':
        log_setlevel(WARNING)
    elif level == 'CRITICAL':
        log_setlevel(CRITICAL)
    elif level == 'RESULT':
        log_setlevel(RESULT)
    elif level == 'ERROR':
        log_setlevel(ERROR)
    else:
        py_printf('Cannot set log level to unsupported log level %s', str(level))


