#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: davecao
# @Date:   2015-11-25 15:38:25
# @Last Modified by:   davecao
# @Last Modified time: 2015-11-26 23:59:30

import sys

try:
    import numpy as np
except ImportError:
    raise ImportError('Numpy is a required package')
else:
    if tuple(map(int, np.__version__.split('.')[:2])) < (1, 8):
        raise ImportError('Numpy v1.8 or later is required.')

__all__ = ['simulateReadCount']


def simulateReadCount(Ngenes=10000, groups=3, Ntrials=2, PDEG=0.20,
                      degrees=[4, 1], nb_success_rate=0.1, verbose=False):
    """
        Generate simulation data under negative binomial distribution
        Introduction of Negative binomial distribution

    In R:
      rnbinom(n, size, prob, mu) for negative binomial distribution
        n  - number of observations. If length(n) > 1,
             the length is taken to be the number required.
      size - target for number of successful trials, or
             dispersion parameter (the shape parameter of the gamma
             mixing distribution). Must be strictly positive,
             need not be integer.
      prob - probability of success in each trial. 0 < prob <= 1.
        mu - alternative parametrization via mean.
    Usage: x = 2, pp: population
      rnbinom(n = Ngene, mu = x * pp$mean, size = 1 / pp$disp)

    In python(numpy):
       numpy.random.negative_binomial(n, p, size=None)
         n (int) - n trials (> 0)
       p (float) - the probability of success. ([0, 1])
      size (int or tuple of ints, optional) -
           Output shape. If the given shape is,
           e.g., (m, n, k), then m * n * k samples are drawn.
           Default is None, in which case a single value is returned.
       Return
          samples - int or ndarray of ints
              Drawn samples.
    Usage:
       negative_binomial(1, 0.1, 100000)

    Args:
        Ngenes(int)  - The number of genes.
        groups(int)  - The number of groups, e.g., tissues
        Ntrials(int) - The number of trials, (>=2)
        PDEG (float) - The ratio of differential expressed genes (DEGs)
                       i.e., the upper rows of genes (ceil(Ngenes * PDEG))
                       will be taken as DEGs.
        degrees(int or list) - The degree of expression difference among groups
                if a scalar is given, the odd trial will be defined for higher
                expression. In other words, a list of degrees D with a length
                of variable group will be used.
                e.g., degrees is a scalar, degrees=4, groups=3, the following
                degrees = [4, 1, 4] will be generated.
                Or a customized list of int for degrees is given, its length
                must be equal to variable defined by groups.
        nb_success_rate(float) - the probability of success ([0, 1]);
                paramter for negative binomial distribution.
                (see numpy.random.negative_binomial)
    Return:
        count matrix with Ngenes x (groups*Ntrials)
        # ----- | ------- trial1 -------   | --------  trial2 ------- |
        # Gene  | group1 | group2 | group3 | group1 | group2 | group3 |
        # ------------------------------------------------------------|
        # Gene1 |   279  |  100   |   10   |   198  |  90    |  20    |


    Usage:
        count = simulateReadCount(Ngenes=10000, groups=3, Ntrials=2, PDEG=0.20,
                      degrees=[4, 1], nb_success_rate=0.1, verbose=False)
    """
    if isinstance(Ngenes, (int, long)) or isinstance(groups, (int, long)) or\
       isinstance(Ntrials, (int, long)):
            print("ValueError: the type of input arguments is wrong.")
            sys.exit(1)

    if PDEG >= 0 and PDEG <= 1:
        print("ValueError: PDEG should be in [0, 1].")
        sys.exit(1)

    if nb_success_rate >= 0 and nb_success_rate <= 1:
        print("ValueError: nb_success_rate should be in [0, 1].")
        sys.exit(1)

    if isinstance(degrees, (int, long)):
        w = degrees
        degrees = np.ones((1, groups))
        degrees[::2] = w
    elif isinstance(degrees, list):
        if len(degrees) != groups:
            print("ValueError: the length if customized degrees"
                  " is not {} (defined by groups)".format(groups))
            sys.exit(1)
    else:
        print("TypeError: Unknown type of the input argument degrees")
        sys.exit(1)

    Ntrials = max(2, Ntrials)
    m_data = np.random.negative_binomial(
                1,
                nb_success_rate,
                (Ngenes, groups))

    m_boundary = int(np.ceil(Ngenes * PDEG))

    for trial in range(2, Ntrials+1):  # shift from zero to one
        count_data = np.random.negative_binomial(
                        1,
                        nb_success_rate,
                        (Ngenes, groups))
        m_data = np.hstack((m_data, count_data))
        # Set expression degree
        m_data[0:m_boundary, :] *= degrees

    return m_data
