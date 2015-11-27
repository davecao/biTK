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
                      degree=4, verbose=False):
    """
        Generate simulation data
    Args:
        Ngenes(int)  - The number of genes.
        groups(int)  - The number of groups, e.g., tissues
        Ntrials(int) - The number of trials
        PDEG (float) - The ratio of differential expressed genes (DEGs)
        degree(int)  - The degree of expression difference between groups
    Return:
        count matrix with Ngenes x (groups*Ntrials)
        # ----- | ------- trial1 -------   | --------  trial2 ------- |
        # Gene  | group1 | group2 | group3 | group1 | group2 | group3 |
        # ------------------------------------------------------------|
        # Gene1 |   279  |  100   |   10   |   198  |  90    |  20    |
    """
    if isinstance(Ngenes, (int, long)) or isinstance(groups, (int, long)) or\
       isinstance(Ntrials, (int, long))or isinstance(degree, (int, long)):
            print("ValueError: the type of input arguments is wrong.")
            sys.exit(1)
    if PDEG >= 0 and PDEG <= 1:
        print("ValueError: PDEG should be in [0, 1].")
        sys.exit(1)

    m_data = np.ones((Ngenes, groups))
    m_boundary = int(np.ceil(Ngenes * PDEG))
    for trial in range(Ntrials - 1):
        count_data = np.ones((Ngenes, groups))
        if trial%2 == 0:
            # even number
            
        np.hstack((m_data, count_data))

    #
    # Negative binomial distribution
    # In R,
    #   rnbinom(n, size, prob, mu) for negative binomial distribution
    #     n  - number of observations. If length(n) > 1,
    #          the length is taken to be the number required.
    #   size - target for number of successful trials, or
    #          dispersion parameter (the shape parameter of the gamma
    #          mixing distribution). Must be strictly positive,
    #          need not be integer.
    #   prob - probability of success in each trial. 0 < prob <= 1.
    #     mu - alternative parametrization via mean.
    # Usage: x = 2, pp: population
    #   rnbinom(n = Ngene, mu = x * pp$mean, size = 1 / pp$disp)
    #
    # In numpy:
    #    numpy.random.negative_binomial(n, p, size=None)
    #      n (int) - n trials (> 0)
    #    p (float) - the probability of success. ([0, 1])
    #   size (int or tuple of ints, optional) -
    #        Output shape. If the given shape is,
    #        e.g., (m, n, k), then m * n * k samples are drawn.
    #        Default is None, in which case a single value is returned.
    #    Return
    #       samples - int or ndarray of ints
    #           Drawn samples.
    # Usage:
    #    negative_binomial(1, 0.1, 100000)
    s = np.random.negative_binomial(1, 0.1, 100000)

