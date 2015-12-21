# -*- coding: utf-8 -*-

# functools is taken as submodule
# biTK.ngs.statistics.functools.{mean, median, ...}
# __all__ = ['functools']

__all__ = []

# biTK.ngs.statistics.{mean, median, ... }
# but biTK.ngs.statistics.functools still hold a copy
#
from biTK.ngs.statistics.functools import *

# push functions in functools.py under functools
# i.e., biTK.ngs.statistics.functools
# from . import functools
# from functools import *

from biTK.ngs.statistics.distributions import *

from biTK.ngs.statistics.normalize import *

from biTK.ngs.statistics.sampling import *
