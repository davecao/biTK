# -*- coding: utf-8 -*-
# The following functions are refered to itertools implemented in pure python
#
import sys
from biTK import PY3K

if PY3K:
    #On Python 3, this will be a zip_longest
    from itertools import zip_longest as zl
else:
    from itertools import izip_longest as zl

__all__ = ['chain', 'product', 'repeat', 'permutations',
        'combinations', 'combinations_with_replacement', 
        'combinations_with_full', 'grouper']

def chain(*iterable):
    # generators
    for it in iterable:
        for inx, element in enumerate(it):
            yield (inx, element)

def product(*args, **kwargs):
    # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    pools = map(tuple, args) * kwargs.get('repeat', 1)
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
    for prod in result:
        yield tuple(prod)


def repeat(object, times=None):
    # repeat(10, 3) --> 10 10 10
    if times is None:
        while True:
            yield object
    else:
        for i in xrange(times):
            yield object

def permutations(iterable, r=None):
    pool = tuple(iterable)
    n = len(pool)
    r = n if r is None else r
    for indices in product(range(n), repeat=r):
        if len(set(indices)) == r:
            yield tuple(pool[i] for i in indices)

def combinations(iterable, r):
    """
     ('ACGT', 2) -> 
        - , ('A', 'C'), ('A', 'G'), ('A', 'T'), 
                -     , ('C', 'G'), ('C', 'T')
                            -     , ('G', 'T')
                                         -  
    """
    pool = tuple(iterable)
    n = len(pool)
    for indices in permutations(range(n), r):
        if sorted(indices) == list(indices):
            yield tuple(pool[i] for i in indices)

def combinations_with_replacement(iterable, r):
    """
     ('ACGT', 2) -> 
        ('A', 'A'), ('A', 'C'), ('A', 'G'), ('A', 'T'), 
                    ('C', 'C'), ('C', 'G'), ('C', 'T')
                                ('G', 'G'), ('G', 'T')
                                            ('T', 'T')
    """
    pool = tuple(iterable)
    n = len(pool)
    for indices in product(range(n), repeat=r):
        if sorted(indices) == list(indices):
            yield tuple(pool[i] for i in indices)

def combinations_with_full(iterable, r):
    """
     ('ACGT', 2) -> 
        ('A', 'A'), ('A', 'C'), ('A', 'G'), ('A', 'T'), 
        ('C', 'A'), ('C', 'C'), ('C', 'G'), ('C', 'T')
        ('G', 'A'), ('G', 'C'), ('G', 'G'), ('G', 'T')
        ('T', 'A'), ('T', 'C'), ('T', 'G'), ('T', 'T')
    """
    pool = tuple(iterable)
    n = len(pool)
    for indices in product(range(n), repeat=r):
        yield tuple(pool[i] for i in indices)

def grouper(n, iterable, opts=(), fillvalue=None):
    """
        x = [ 1, 2, 3, 4, 5]
        grouper(2, x) -->[(1,2), (3,4), (5, None)]
    """
    if not isinstance(opts, tuple):
        print("opts in grouper should be a tuple")
        sys.exit(1)
    args = [iter(iterable)] * n
    #return zl(fillvalue=fillvalue, *args)
    for x in zl(fillvalue=fillvalue, *args):
        r = x + opts
        #print(r)
        #sys.stdout.flush()
        yield r
