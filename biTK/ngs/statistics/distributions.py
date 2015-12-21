#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: davecao
# @Date:   2015-11-25 15:38:25
# @Last Modified by:   davecao
# @Last Modified time: 2015-12-21 15:42:55

import os
import sys
import biTK
from biTK.ngs.statistics.normalize import normalize_by_median
from biTK.ngs.utils import product
from biTK import PY3K

if PY3K:
    # On Python 3, this will be a zip_longest
    from itertools import zip_longest as zl
else:
    from itertools import izip_longest as zl

try:
    import numpy as np
except ImportError:
    raise ImportError('Numpy is a required package')
else:
    if tuple(map(int, np.__version__.split('.')[:2])) < (1, 8):
        raise ImportError('Numpy v1.8 or later is required.')

__all__ = ['simulateReadCount', 'negative_binomial',
           'getExperimentData']

arab_file = biTK.DATA + os.sep + "arab.txt"


def getExperimentData(filename,
                      skip_header=1,
                      usecols=(1, 2, 3, 4, 5, 6), **kwargs):
    """ Load arab.txt """
    # colname = ["gene_id", "mock1", "mock2", "mock3",
    #           "hrcc1", "hrcc2", "hrcc3"]
    replicates = kwargs.get('replicates', 3)
    delimiter = kwargs.get('delimiter', ' ')
    ax = kwargs.get('ax', 1)
    data = np.genfromtxt(filename, skip_header=1,
                         delimiter=delimiter, usecols=usecols)
    # split into groups
    nrows, ncols = data.shape
    # normalize the data by columns: this will eliminate errors
    norm_data = normalize_by_median(data, replicates=replicates)
    # split by replicates, keep rows
    dim = ncols/replicates
    r_data = np.hsplit(norm_data, dim)
    # r_data = np.asarray(r_data)
    # mean, var of subarrays
    result_mu = []
    result_r = []
    for i in range(dim):  # do in-place
        #  E = mu;
        #  Var = mu + phi*mu^{2}
        #  phi = 1/r
        #
        #  prob:(success)
        #            var - E           E^{2}
        #       p = ---------- , r = -----------
        #              Var            Var - E
        # Remove the rows which only contains zeros
        inx = ~np.all(r_data[i] == 0, axis=ax)
        m = np.mean(r_data[i][inx], axis=ax)
        v = np.var(r_data[i][inx], axis=ax)
        # p = (v-m)/v
        r = m*m / (v - m)
        # remove negatives in r
        r_inx = r > 0
        # save mean and r
        result_mu.extend(m[r_inx])
        result_r.extend(r[r_inx])
    # result: 2xn array
    #   0xn: [mu, mean]
    #   1xn: [r=1/phi]
    result = np.vstack((result_mu, result_r))
    return result


def poisson_gamma_mixture(shape, scale=1.0, size=None):
    """
        Generate an integer from Poisson gamma mixture distribution

    Args:
        shape (scalar>0): The shape of the gamma distribution.
                          "k" is used sometimes.
        scale (scalar>0, optional): the scale of the gamma distribution
                         "theta" is used sometimes.
        see numpy.random.gamma
    Returns:
        Draw a sample from Poisson gamma mixture distribution
    """
    lam = np.random.gamma(shape, scale=scale, size=size)
    s = np.random.poisson(lam=lam, size=size)
    return s


def negative_binomial(mu, r, size=None):
    """
        Generate a sample from negative binomial distribution

    Args:
        mu (float or array): mean value, if mu or r are two arrays of
                             same length, random variables will be generated
                             by each pair of mu and r.
        r(float or array): success trials.
        size (integer or tuple): the number of samples to generate.
                    if size is a tuple of two elements, (nrows, ncols)
                     (mu, r)     (nrows, ncols)
                     0.5, 1 -->   ncols from negative binomial(mu,r)
                     0.5, 2 -->
                     0.5, 4 -->
                if nrows > size[0], enable to choice with repeated pair.
                  using np.random.choice(size, nrows, repeat=True)
                else nrows < size[0], disable to choice with repeated pair.
                  using np.random.choice(size, nrows, repeat=False)
        e.g.,
                mu = 0.5 and r = [1, 2, 4]
            The combination of mu and r will be used as followed:
                mu = [0.5, 0.5, 0.5]
                r = [1, 2, 4]
            if size (integer) is given, then repeat the times of size.

    Returns:
        A n-elements array of integers
    """
    if size is not None:
        if not (isinstance(size, int) or isinstance(size, tuple)):
            print("ValueError: size should be an integer or a tuple")
            sys.exit(1)

    params_ = []

    if isinstance(mu, float) and isinstance(r, float):
        return poisson_gamma_mixture(r, scale=mu/r, size=size)
    elif isinstance(mu, float) and isinstance(r, np.ndarray):
        for m_val, r_val in product([mu], r):
            params_.append([m_val, r_val])
    elif isinstance(mu, np.ndarray) and isinstance(r, float):
        for m_val, r_val in product(mu, [r]):
            params_.append([m_val, r_val])

#        if size is None:
#            s = [poisson_gamma_mixture(r_val, scale=m_val/r_val)
#                 for m_val, r_val in product(mu, [r])]
#        else:
#            s = [[poisson_gamma_mixture(r_val, scale=m_val/r_val)
#                 for m_val, r_val in product(mu, [r])] for i in range(size)]
#        return s
    elif isinstance(mu, np.ndarray) and isinstance(r, np.ndarray):
        if mu.shape != r.shape:
            print("ValueError: the shape of mu and r should be same.")
            sys.exit(1)
        for m_val, r_val in zl(mu, r):
            params_.append([m_val, r_val])

#        if size is None:
#            s = [poisson_gamma_mixture(r_val, scale=m_val/r_val)
#                 for m_val, r_val in zl(mu, r)]
#        else:
#            s = [[poisson_gamma_mixture(r_val, scale=m_val/r_val)
#                 for m_val, r_val in zl(mu, r)] for i in range(size)]
#        return s
    else:
        print("ValueError: Wrong types of input arguments")
        sys.exit(1)
    # convert to numpy array
    params_ = np.asarray(params_)
    if size is None:
        # draw samples, the number of which is the number of pairs of mu and r
        return [poisson_gamma_mixture(r_val, scale=m_val/r_val)
                for m_val, r_val in params_]
    else:
        # preallocate
        s = np.zeros(size)
        # randomly select pairs
        rep = False
        params_len = len(params_)
        # nrows, ncols = size
        if params_len < size[0]:
            rep = True
        p_inx = np.random.choice(params_len, size[0], replace=rep)
        params_pair = params_[p_inx]
        if s.ndim == 1:
            s[:] = [poisson_gamma_mixture(r_val, scale=m_val/r_val)
                    for m_val, r_val in params_pair]
        elif s.ndim == 2:
            s[:, :] = [[
                poisson_gamma_mixture(
                    r_val, scale=m_val/r_val, size=s.shape[1])
                for m_val, r_val in params_pair] for i in range(s.shape[0])]
        else:
            print("Only 2d array supported at present. Now return a zeros")
        return s


def simulateReadCount(Ngenes=10000, groups=3, Ntrials=3, PDEG=0.20,
                      degrees=[4, 1, 1], nb_success_rate=0.1, verbose=False):
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
        # ----- | ------- time 1 -------   | --------  time 2 ------- |
        # Gene  |  exp1  |  exp2  |   exp3 |   exp1 |   exp2 |  exp3  |
        # ------------------------------------------------------------|
        # Gene1 |   279  |  100   |   10   |   198  |  90    |  20    |


    Usage:
        count = simulateReadCount(Ngenes=10000, groups=3, Ntrials=2, PDEG=0.20,
                      degrees=[4, 1], nb_success_rate=0.1, verbose=False)
        On terminal
        python -c "import sys;import biTK;import numpy as np; c=biTK.ngs.statistics.simulateReadCount();np.savetxt(sys.stdout, c, fmt='%5d')"
    """
    if not isinstance(Ngenes, (int, long)) or \
       not isinstance(groups, (int, long)) or \
       not isinstance(Ntrials, (int, long)):
            print("ValueError: the type of input arguments is wrong.")
            sys.exit(1)

    if not (PDEG >= 0 and PDEG <= 1):
        print("ValueError: PDEG should be in [0, 1].")
        sys.exit(1)

    if not (nb_success_rate >= 0 and nb_success_rate <= 1):
        print("ValueError: nb_success_rate should be in [0, 1].")
        sys.exit(1)

    if isinstance(degrees, (int, long)):
        w = degrees
        degrees = np.ones((1, groups))
        degrees[::2] = w
    elif isinstance(degrees, list):
        degrees = np.asarray(degrees)
        if len(degrees) != groups:
            print("ValueError: the length if customized degrees"
                  " is not {} (defined by groups)".format(groups))
            sys.exit(1)
    else:
        print("TypeError: Unknown type of the input argument degrees")
        sys.exit(1)

    # Generate simulate data
    # 1. load arab data
    mu, r = getExperimentData(arab_file)
    Ntrials = max(2, Ntrials)
    # 2. allocate a matrix
    #   Note: ncols = groups * Ntrials
    #         nrows = Ngenes

    # m_data = np.random.negative_binomial(
    #            1,
    #            nb_success_rate,
    #            (Ngenes, groups))
    #
    m_data = negative_binomial(mu, r, size=(Ngenes, groups*Ntrials))
    m_boundary = int(np.ceil(Ngenes * PDEG))

    for trial in range(2, Ntrials+1):  # shift from zero to one
        count_data = np.random.negative_binomial(
                        1,
                        nb_success_rate,
                        (Ngenes, groups))
        # Set expression degree
        count_data[0:m_boundary, :] *= degrees
        # Combine data
        m_data = np.hstack((m_data, count_data))

    return m_data
