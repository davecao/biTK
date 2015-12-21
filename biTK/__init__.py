# -*- coding: utf-8 -*-

__author__ = "Wei Cao"
__contact__ = "davecao@bi.a.u-tokyo.ac.jp"
__date__ = "2015/09/11"
__version__ = '1.0.0'

import os
import os.path
import sys
import warnings
import sysconfig
import types

if sys.version_info[:2] < (2, 7):
    raise Exception('Python version less than 2.7')
VERSION = __version__
#try:
#    import numpy as np
#except ImportError:
#    raise ImportError('Numpy is a required package')
#else:
#    if tuple(map(int, np.__version__.split('.')[:2])) < (1, 8):
#        raise ImportError('Numpy v1.8 or later is required.')

from jinja2 import Environment, FileSystemLoader

_PY3K = PY3K = sys.version_info[0] > 2
PY2K = not PY3K


def is_float(n):
    return isinstance(n, float)

TPLPATH = "{}/ngs/templates/".format(os.path.dirname(__file__))
jinja2_ENV = Environment(loader=FileSystemLoader(TPLPATH))
jinja2_ENV.add_extension('jinja2.ext.loopcontrols')
jinja2_ENV.tests['Float'] = is_float
DATA = "{}/ngs/data".format(os.path.dirname(__file__))

__all__ = []

if PY3K:
    from copyreg import __newobj__ as reduce_newobj
    string_types = str,
    integer_types = int,
    class_types = type,
    text_type = str,
    binary_type = bytes
else:
    from copy_reg import __newobj__ as reduce_newobj
    string_types = basestring,
    integer_types = (int, long),
    class_types = (type, types.ClassType),
    text_type = unicode,
    binary_type = str

from biTK.ngs import *

