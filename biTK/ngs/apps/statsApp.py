#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import uuid
import time
import copy

from signal import signal, SIGPIPE, SIG_DFL
from optparse import OptionParser, OptionGroup
from collections import Counter
from operator import add
#from functools import wraps

from biTK.ngs.sequence import nucleic_alphabet
from biTK.ngs.io import AlignIO
from biTK.ngs.io.pathtools import isReadable, isfile, isdir
from biTK.ngs.statistics import mean, median, quantile
from biTK.ngs.utils import Console, combinations_with_full, \
                            AnalysisItems, IndentedHelpFormatterWithNL, \
                            StatsAnalysisReport
from biTK import string_types
#from multiprocessing.pool import ThreadPool as Pool
from biTK.ngs.concurrent import ThreadPool
from multiprocessing import cpu_count

try:
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import *
    matplotlib_version = matplotlib.__version__
except ImportError:
    raise ImportError('Matplotlib is a required package')


__author__ = "Wei Cao"
__contact__ = "davecao@bi.a.u-tokyo.ac.jp"
__date__ = "2015/08/27"
__version__ = "0.1"
__copyright__ = """
    Feel free to do whatever you like with this code.
    """

# reset
signal(SIGPIPE, SIG_DFL)
# global
LOGGER = None
ThPool = None
#ACTIONS = {
#    'stats': statistics_report, 
#    'tbs_head': tbs_head, 
#    'tbs_tail': tbs_tail, 
#    'tbl_head': tbl_head,
#    'tbl_tail': tbl_tail,
#    'save'    : save2fastq
#}
ACTIONS = ['stats', 'tbs_head', 'tbs_tail', 'tbl_head', 'tbl_tail', 'save']
## Cache decorator
class memoize:
    # from http://avinashv.net/2008/04/python-decorators-syntactic-sugar/
    def __init__(self, function):
        self.function = function
        self.memoized = {}

    def __call__(self, *args):
        try:
            return self.memoized[args]
        except KeyError:
            self.memoized[args] = self.function(*args)
            return self.memoized[args]

def load_data(filename, minLen, 
                io_plug = 'FastQIO_multithread',
                fastq_type = 'single',
                compressed=None, 
                alphabet=nucleic_alphabet, 
                quality_score_fmt='phred33',
                nthreads = None,
                chks = None,
                verbose=False):
    """
        Read NGS data in fastq format 
    """
    if nthreads is None:
        nthreads = cpu_count()
    parser = AlignIO(filename, ConcreteIO=io_plug, compressed=compressed)
    seqslist_full = parser.parse(alphabet=alphabet, 
                                 quality_score_fmt=quality_score_fmt,
                                 nthreads=nthreads, chunksize=chks, 
                                 verbose=verbose)
    seqslist = []
    removed_list = []
    s_ap = seqslist.append
    rl_ap = removed_list.append

    max_length = 0
    if minLen:
        for seq in seqslist_full:
            length = seq.length()
            max_length = max(length, max_length)
            if length >= minLen:
                #seqslist.append(seq)
                s_ap(seq)
            else:
                #removed_list.append(seq)
                rl_ap(seq)
    else:
        for seq in seqslist_full:
            length = seq.length()
            max_length = max(length, max_length)
        seqslist = seqslist_full
    return seqslist, removed_list, max_length

def gc_content_all(contents, alphabet=nucleic_alphabet, pth="."):
    """
        Investigate GC content of input sequences
    """
    #contents = [0] * len(alphabet)
    #for seq in seqslist:
    #    seq_cons = seq.tally(alphabet=alphabet)
    #    contents = map(add, seq_cons, contents)
    #nt = dict(zip(list(alphabet), contents))
    #GC_content = float(nt['G'] + nt['C'])/sum(contents)*100.0

    from itertools import chain
    data = list(chain.from_iterable(contents))
    CG_inx = [alphabet.index('C'), alphabet.index('G')]
    GC_sum = 0.0
    for item in chain(contents):
        selected_ = [item[i] for i in CG_inx]
        GC_sum += sum(selected_, 0.0)

    GC_content = GC_sum/sum(data, 0.0) * 100
    exp = AnalysisItems(
            name="GC contents", 
            description="GC contents in total: {:.3f}%".format(GC_content), 
            data=GC_content, 
            hasTable=False, 
            hasImg=False,
            img_path=None, 
            img_name=None)
    return exp, GC_content

def boxWisherPlot(data, mean_data, median_data, quantile_data, fliersMarker='',
                figsize=(10,6),
                title=None, xlabel=None, ylabel=None, xygrid=[False, True],
                filledColors = 'yellow',
                notch=0, sym='+', vert=1, whis=1.5, 
                *args, **kwargs):
    """
       customize matplotlib's boxplot
        see http://matplotlib.org/examples/pylab_examples/boxplot_demo2.html
    
    Args:
        data (list of list): 
        mean_data (list):
        median_data(list):
        quantile_data(list): 
        fliersMarker (str) : default is '+'
        figsize (tuple):  set figure size, default is (10,6)
        title (str) : the title of the figure.
        xlabel (str): set x axis label 
        ylabel (str): set y axis label
        xygrid (list) : default is [False, True], corresponding to enable or 
                        disable x-grid and y-grid.
        filledColors (str or list): the color used to fill the boxes 
        notch (bool) : default is 0 
        sym (str) :  default is '+' 
        vert : default is 1 
        whis : default is 1.5
        other see matplotlib's boxplot()
    """
    from itertools import chain

    n_box = len(data)
    # check dimensions
    ndims = len(set([n_box, len(mean_data),len(median_data), len(quantile_data)]))
    if ndims != 1:
        print("Input: data, mean, median, quantile incompatible")
        sys.exit(1)

    fig, ax = plt.subplots(figsize=figsize)
    plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

    bp = plt.boxplot(data, notch=notch, sym=sym, vert=vert, whis=whis, 
                    *args, **kwargs)

    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='red', marker='+')

    # Add grids
    if xygrid[0]:
        # vertical grids
        ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
              alpha=0.5)
    if xygrid[1]:
        # horizontal grids
        ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
              alpha=0.5)
    # Hide these grid behind plot objects
    ax.set_axisbelow(True)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # Now fill the boxes with desired colors
    # e.g. ['darkkhaki','royalblue'] 
    boxColors = [0]*n_box
    if isinstance(filledColors, string_types):
        boxColors = [ filledColors for i in boxColors ]
    elif isinstance(filledColors, list):
        counter_f=0
        n_filledColors = len(filledColors)
        for i in range(n_box):
            if counter_f < n_filledColors:
                boxColors[i] = filledColors[counter_f]
                counter_f += 1
            else:
                counter_f = 0
                boxColors[i] = filledColors[counter_f]

    # modify the boxes
    for i in range(n_box):
        # quantile data
        (q1, q2, q3, q4) = quantile_data[i]
        # lower cap
        bp['caps'][i*2].set_ydata([q1, q1])
        # higher cap
        bp['caps'][i*2 + 1].set_ydata((q4, q4))
        # lower whisker
        bp['whiskers'][i*2].set_ydata([q1, q2])
        # higher whisker
        bp['whiskers'][i*2 + 1].set_ydata((q3, q4))

        box = bp['boxes'][i]
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        boxCoords = zip(boxX, boxY)
        #Filled with colors
        boxPolygon = Polygon(boxCoords, facecolor=boxColors[i])
        ax.add_patch(boxPolygon)
        # Now draw the median lines back over what we just filled in
        med = bp['medians'][i]
        #medianX = []
        #medianY = []
        #for j in range(2):
        #    medianX.append(med.get_xdata()[j])
            #medianY.append(med.get_ydata()[j])
        #    plt.plot(medianX, median_data, 'k')
        #medians[i] = medianY[0]
        # Finally, overplot the sample averages, with horizontal alignment
        # in the center of each box
        plt.plot(mean(med.get_xdata()), mean_data[i],
                    color='w', marker='*', markeredgecolor='k')
    # Set the axes ranges and axes labels
    ax.set_xlim(0.5, n_box+0.5)
    top = 40
    bottom = -5
    ax.set_ylim(bottom, top)
    xtickNames = plt.setp(ax, xticklabels=range(1,n_box+1))
    #plt.setp(xtickNames, rotation=45, fontsize=8)
    plt.setp(xtickNames, fontsize=8)

    return fig

def sequence_quality_per_base(scores_mat, img_file=None, dpi=300, pth ="."):
    """
        Generate data for BoxWhisker plot 
    """
    name = u"Per base sequence quality"
    description = u"An overview of the range of quality values across all bases at each position"
    if img_file is None:
        img_file = 'sequence_quality_per_base.png'
        img_path = pth + os.sep + img_file
    else:
        img_path = pth + os.sep + img_file
    #from itertools import chain
    #scores_mat = []
    #for seq in seqslist:
    #    scores_mat.append(seq.scores)
    # row is position in read
    # col is scores
    scores_mat_trans = map(list, zip(*scores_mat))
    n_pos = len(scores_mat_trans)
    result = [0] * n_pos
    mean_all = [0] * n_pos
    median_all = [0] * n_pos
    fliers_xy = [0] * n_pos

    for inx, item in enumerate(scores_mat_trans):
        item_sorted = sorted(item)
        mean_val = mean(item_sorted)
        median_val = median(item_sorted)
        quantile_val = quantile(item_sorted, prob=[0.1, 0.25, 0.75, 0.9], 
                                alphap=.4, betap=.4)
        mean_all[inx] = mean_val
        median_all[inx] = median_val
        result[inx] = quantile_val
 
    fig = plt.figure()
    fig.subplots_adjust(wspace=0.22,hspace=0.2)
    ax = fig.add_subplot(111)
    fig = boxWisherPlot(scores_mat_trans, mean_all, median_all, result,
            title='Per base sequence quality',
            xlabel='Position in read',
            ylabel='Scores',
            xygrid=[False, True])

    fig.savefig(img_path, dpi=dpi) # save to file
    plt.close(fig) # close the figure
    exp = AnalysisItems(
            name=name, 
            description=description, 
            data=scores_mat_trans, 
            hasTable=False, 
            hasImg=True,
            img_path=img_path, 
            img_name=img_file)
    return exp

def quality_scores_per_sequence(quality_scores_mat, img_file=None, dpi=300, pth = "."):
    """
        Mean sequence quality(Phred score) vs occurencies
    Args:
        quality_scores_mat (list of list):
            row: number of sequences
            col: scores on each position
    """
    name = u"Per sequence quality scores"
    description = u"Mean sequence quality(Phred score) vs occurencies"
    if img_file is None:
        img_file = 'quality_scores_per_sequence.png'
        img_path = pth + os.sep + img_file
    else:
        img_path = pth + os.sep + img_file

    from itertools import chain
    #data = list(chain.from_iterable(quality_scores_mat))
    data = []
    for item in chain(quality_scores_mat):
        data.append(sum(item, 0.0) / len(item))
    c = Counter(data)
    c_sorted = sorted(zip(c.keys(), c.values()))
    
    fig = plt.figure()
    fig.subplots_adjust(wspace=0.22,hspace=0.2)
    ax = fig.add_subplot(111)
    ax.plot(*zip(*c_sorted))
    ax.set_title('Per sequence quality scores')
    ax.set_xlabel('Mean sequence quality')
    ax.set_ylabel('Number of reads')
    ax.yaxis.grid(True)
    fig.savefig(img_path, dpi=dpi) # save to file
    plt.close(fig) # close the figure
    exp = AnalysisItems(
            name=name, 
            description=description, 
            data=c_sorted, 
            hasTable=False, hasImg=True,
            img_path=img_path, 
            img_name=img_file)
    return exp

def sequence_content_per_base(contents, alphabet=nucleic_alphabet, img_file=None, dpi=300, pth = "."):
    """
        contents(%) of A,T, G and C on each position
    Args:
        contents (list of list):
            row: each position in read (bp)
            col: nucleic_alphabet order
    """
    name = u"Per base sequence contents"
    description = u"contents(%) of A,T, G and C on each positioin"
    if img_file is None:
        img_file = 'sequence_content_per_base.png'
        img_path = pth + os.sep + img_file
    else:
        img_path = pth + os.sep + img_file

    ACGT_inx = [alphabet.index('A'), alphabet.index('C'),
                alphabet.index('G'), alphabet.index('T')]
    from itertools import chain
    c = []
    for item in chain(contents):
        selected_ = [item[i] for i in ACGT_inx]
        #s = sum(selected_, 0.0)
        s = sum(item, 0.0)
        if s == 0:
            s = 1
        c.append([d/s*100 for d in selected_])

    fig = plt.figure()
    fig.subplots_adjust(wspace=0.22,hspace=0.2)
    ax = fig.add_subplot(111)
    c_trans = map(list, zip(*c))
    for item in chain(c_trans):
        ax.plot(item)
    ax.set_title('Per base sequence contents')
    ax.set_xlabel('Position in read (bp)')
    ax.set_ylabel('contents(%)')
    ax.legend(['A', 'C', 'G', 'T'],loc='best', fancybox=True, 
                shadow=False, framealpha=0.5)
    ax.yaxis.grid(True)
    fig.savefig(img_path, dpi=dpi) # save to file
    plt.close(fig) # close the figure
    exp = AnalysisItems(
            name=name, 
            description=description, 
            data=c_trans, 
            hasTable=False, hasImg=True,
            img_path=img_path, 
            img_name=img_file)
    return exp

def gc_content_per_base(contents, alphabet=nucleic_alphabet, img_file=None, dpi=300, pth="."):
    """
        GC content (%) on each position
    Args:
        contents (list of list) :
            row: each position in read (bp)
            col: nucleic_alphabet order
    """
    name = u"Per base gc contents"
    description = u"GC contents(%) on each position"
    if img_file is None:
        img_file = 'gc_content_per_base.png'
        img_path = pth + os.sep + img_file
    else:
        img_path = pth + os.sep + img_file
    GC_inx = [alphabet.index('G'),alphabet.index('C')]
    from itertools import chain
    gc = []
    for item in chain(contents):
        selected_ = sum([item[i] for i in GC_inx])
        s = sum(item, 0.0)
        if s == 0:
            s = 1
        gc_row = selected_/s*100
        gc.append(gc_row)
    fig = plt.figure()
    fig.subplots_adjust(wspace=0.22,hspace=0.2)
    ax = fig.add_subplot(111)
    ax.plot(gc)
    ax.set_title('Per base GC contents')
    ax.set_xlabel('Position in read (bp)')
    ax.set_ylabel('GC contents(%)')
    ax.yaxis.grid(True)
    fig.savefig(img_path, dpi=dpi) # save to file
    plt.close(fig) # close the figure
    exp = AnalysisItems(
            name=name, 
            description=description, 
            data=gc, hasTable=False, hasImg=True,
            img_path=img_path, 
            img_name=img_file)
    return exp

def gc_content_per_sequence(seqslist, alphabet=nucleic_alphabet, img_file=None, dpi=300, pth="."):
    """
        distribution of GC content (%) per sequence
    Args:
        seqslist (list of Sequence)
    """
    name = u"Per sequence gc contents"
    description = u"GC contents(%) distribution from each sequence"
    if img_file is None:
        img_file = 'gc_content_per_sequence.png'
        img_path = pth + os.sep + img_file
    else:
        img_path = pth + os.sep + img_file
    gc = [0]*len(seqslist)
    for inx, seq in enumerate(seqslist):
        gc[inx] = seq.GC_content
        #gc.append(seq.get_GC_content(alphabet=alphabet))
    gc_count = Counter(gc)
    gc_sorted = sorted(zip(gc_count.keys(), gc_count.values()))

    fig = plt.figure()
    fig.subplots_adjust(wspace=0.22,hspace=0.2)
    ax = fig.add_subplot(111)
    ax.plot(*zip(*gc_sorted))
    ax.set_title('Per sequence GC contents')
    ax.set_xlabel('GC contents(%)')
    ax.set_ylabel('Number of reads')
    ax.yaxis.grid(True)
    fig.savefig(img_path, dpi=dpi) # save to file
    plt.close(fig) # close the figure
    exp = AnalysisItems(
            name=name, 
            description=description, 
            data=gc_sorted, hasTable=False, hasImg=True,
            img_path=img_path, 
            img_name=img_file)
    return exp

def N_content_per_base(contents, alphabet=nucleic_alphabet, img_file=None, dpi=300, pth="."):
    """
        N(%) on each position in read(bp)
    """
    name = u"N(%) on each position"
    description = u"N contents(%) on each position"
    if img_file is None:
        img_file = 'N_content_per_base.png'
        img_path = pth + os.sep + img_file
    else:
        img_path = pth + os.sep + img_file
    N_inx = alphabet.index('N')
    from itertools import chain
    N_content = []
    for item in chain(contents):
        d = item[N_inx]
        s = sum(item, 0.0)
        if s == 0:
            s = 1
        N_content.append(d/s*100)

    fig = plt.figure()
    fig.subplots_adjust(wspace=0.22, hspace=0.2)
    ax = fig.add_subplot(111)
    ax.plot(N_content)
    ax.set_title('Per base sequence contents')
    ax.set_xlabel('Position in read (bp)')
    ax.set_ylabel('contents(%)')
    ax.yaxis.grid(True)
    fig.savefig(img_path, dpi=dpi) # save to file
    plt.close(fig) # close the figure
    exp = AnalysisItems(
            name=name, 
            description=description, 
            data=N_content, hasTable=False, hasImg=True,
            img_path=img_path, 
            img_name=img_file)
    return exp

def seq_len_dist(data, img_file=None, dpi=300, pth="."):
    """
        distribution of sequence length
    Args:
        data (list): sequence length
    """
    name = u"Distribution of sequence lengths"
    description = u""
    if img_file is None:
        img_file = 'seq_len_dist.png'
        img_path = pth + os.sep + img_file
    else:
        img_path = pth + os.sep + img_file
    l = Counter(data)
    l_sorted = sorted(zip(l.keys(), l.values()))

    fig = plt.figure()
    fig.subplots_adjust(wspace=0.22,hspace=0.2)
    ax = fig.add_subplot(111)
    ax.plot(*zip(*l_sorted))
    ax.set_title('Distribution of sequence lengths')
    ax.set_xlabel('Sequence length')
    ax.set_ylabel('Number of reads')
    ax.yaxis.grid(True)
    fig.savefig(img_path, dpi=dpi) # save to file
    plt.close(fig) # close the figure
    exp = AnalysisItems(
            name=name, 
            description=description, 
            data=l_sorted, hasTable=False, hasImg=True,
            img_path=img_path, 
            img_name=img_file)
    return exp

def sequence_duplicates(seqslist, pth="."):
    """
        Count duplicate reads
    """
    name = u"Occurrencies of duplicate reads"
    description = u"Count the duplicates"
    L = len(seqslist)
    LOGGER.timeit(label='Counter')
    u_seq = Counter(seqslist)
    LOGGER.report(msg='Counter - Completed in %.2fs.', label='Counter')
    sum_val = sum(u_seq.values())
    # most expensive part: 282.15s
    LOGGER.timeit(label='sort Counter')
    u_sorted = sorted(zip(u_seq.values(),u_seq.keys())) # in increasing order
    u_sorted = [(b, a, float(a)/sum_val*100) for a,b in u_sorted]
    LOGGER.report(msg='sort Counter - Completed in %.2fs.', label='sort Counter')
    exp = AnalysisItems(
            name=name, 
            description=description, 
            data=u_sorted, hasTable=True, tableColName=['Read', 'Occurrence', 'Percentage'], hasImg=False,
            img_path=None, 
            img_name=None)
    return exp
    #return u_seq, u_sorted[-1]

def kmer_nheads_reads(seqslist, k, n, alphabet=None, verbose=False, pth="."):
    """
        Statistics on reads
    Args:
        k (int): k-mer in the reads
        n (int): n bases in the head of reads
    """
    name=u"{}-mer analysis of reads".format(k)
    description=u"{}-mer occurred in all reads".format(k)
    # words count
    n_elems_head = []
    kmer_nt = {}
    # count occurencies of k-mer 
    for seq in seqslist:
        n_elems_head.append(seq.get_n_head(n))
        for mer in seq.words(k):
            try:
                kmer_nt[mer] += 1
            except KeyError:
                kmer_nt[mer] = 1
                #print('KeyError:{}'.format(mer))
    sum_kmer = sum(kmer_nt.values(), 0.0)
    kmer_nt_sorted = sorted(zip(kmer_nt.values(), kmer_nt.keys()), reverse=True)
    kmer_nt_sorted = [ (b, a, a/sum_kmer*100) for a, b in kmer_nt_sorted ]
    # count occurrencies of n bases in the head
    n_heads = Counter(n_elems_head)
    sum_heads = sum(n_heads.values(), 0.0)
    n_heads_sorted = sorted(zip(n_heads.values(), n_heads.keys()), reverse=True)
    n_heads_sorted = [ (b, a, a/sum_heads*100) for a, b in n_heads_sorted ]
    exp = AnalysisItems(
            name=name, 
            description=description, 
            data=kmer_nt_sorted, tableColName=['{}-mer'.format(k), 'Occurrence', 'Percentage'], hasTable=True, hasImg=False,
            img_path=None, 
            img_name=None)
    return exp
    #return kmer_nt_sorted, n_heads

def stats_preprocess(seqslist, npos, alphabet=nucleic_alphabet):
    """
        Count occurencies on each position
    Args:
        seqslist (list of Sequence objs): list of objects of Sequence
        npos (int) : maximum sequence length in seqslist
    """
    nseqs = len(seqslist)
    ncols = len(alphabet)
    length_dist = []
    quality_scores_mat = []
    ld_ap = length_dist.append
    qs_ap = quality_scores_mat.append
    # initial list of list
    # npos: 
    # ncols: alphabet order
    contents = [[ 0 for r in range(ncols) ] for c in range(npos)]

    #for seq in seqslist:
    for i in range(nseqs):
        s = list(seqslist[i].get_raw_seq())
        ld_ap(seqslist[i].length())
        # score is missing, i.e., sequence length is shorter than 
        # the length of the longest sequence
        #scores = [ 0 for i in range(npos) ]
        #for i, v in enumerate(seq.scores):
        #    scores[i] = v
        scores = [0] * npos
        n = len(seqslist[i].scores)
        scores[:n] = seqslist[i].scores
        qs_ap(scores)

        for row, val in enumerate(s):
            col = alphabet.index(val[0])
            contents[row][col] += 1

    return contents, quality_scores_mat, length_dist

def trim_adapter_sequence(adapter):
    """
        Trim adapter sequences from reads
    """
    pass

def tbs_head(seqslist, score, inPlace=False, verbose=False):
    """
        Trim ngs reads by quality scores from heads
    Args:
        seqslist (list of Sequence) : list of Sequence objects.
        score (float)  : The cutoff score for triming sequences
        inPlace (bool) : Trim each sequence inplace if True, default is False.
                         If False, duplicate the list of objects and return a 
                         new list of objects that each sequence was trimmed.
        verbose (bool) : show verbosity
    """
    if inPlace:
        # shallow copy - reference
        seqslist_copy = seqslist
    else:
        # deep copy
        seqslist_copy = copy.deepcopy(seqslist, memo=None)
    for seq in seqslist_copy:
        seq.trim_score_head(score)
    if not inPlace:
        return seqslist_copy

def tbs_tail(seqslist, score, inPlace=False, verbose=False):
    """
        Trim ngs reads by quality scores from tails
    Args:
        seqslist (list of Sequence) : list of Sequence objects.
        score (float)  : The cutoff score for triming sequences
        inPlace (bool) : Trim each sequence inplace if True, default is False.
                         If False, duplicate the list of objects and return a 
                         new list of objects that each sequence was trimmed.
        verbose (bool) : show verbosity
    """
    if inPlace:
        # shallow copy - reference
        seqslist_copy = seqslist
    else:
        # deep copy
        seqslist_copy = copy.deepcopy(seqslist, memo=None)
    for seq in seqslist_copy:
        seq.trim_score_tail(score)
    if not inPlace:
        return seqslist_copy

def tbl_head(seqslist, n, inPlace=False, verbose=False):
    """
        Trim ngs reads by a given length
    Args:
        seqslist (list of Sequence) : list of Sequence objects.
        n (int)  : Trim the n bases in head of sequences
        inPlace (bool) : Trim each sequence inplace if True, default is False.
                         If False, duplicate the list of objects and return a 
                         new list of objects that each sequence was trimmed.
        verbose (bool) : show verbosity
    """
    if inPlace:
        # shallow copy - reference
        seqslist_copy = seqslist
    else:
        # deep copy
        seqslist_copy = copy.deepcopy(seqslist, memo=None)
    for seq in seqslist_copy:
        seq.trim_len_head(n)

    if not inPlace:
        return seqslist_copy

def tbl_tail(seqslist, n, inPlace=False, verbose=False):
    """
        Trim ngs reads by a given length
    Args:
        seqslist (list of Sequence) : list of Sequence objects.
        ntails (int)  : Trim the n bases in tail of sequences
        inPlace (bool) : Trim each sequence inplace if True, default is False.
                         If False, duplicate the list of objects and return a 
                         new list of objects that each sequence was trimmed.
        verbose (bool) : show verbosity
    """
    if inPlace:
        # shallow copy - reference
        seqslist_copy = seqslist
    else:
        # deep copy
        seqslist_copy = copy.deepcopy(seqslist, memo=None)
    for seq in seqslist_copy:
        seq.trim_len_tail(n)
    if not inPlace:
        return seqslist_copy

def statistics_report(seqslist, removed_list, max_length, cmd, opt):
    """
        do statistical analysis
    """
    # Statistics
    inputfile = opt.fastqfile
    input_bsnames = ','.join([os.path.basename(p) for p in inputfile])
    reportdir = opt.reportdir
    n_seqs = len(seqslist)
    n_removed = len(removed_list)

    LOGGER.info("Fastq File: {}".format(inputfile))
    LOGGER.info("Total sequences: {}".format(len(seqslist)+len(removed_list)))
    LOGGER.info("Sequences: {}".format(len(seqslist)))
    trim_min_len = opt.trim_min_len if opt.trim_min_len else 0
    LOGGER.info("Removed sequences(<{}): {}".format(trim_min_len, 
                                                    len(removed_list)))
    LOGGER.info("Nucleic types : {}".format(nucleic_alphabet))

    #stats_ngs_reads(seqslist, opt.stats_kmer, opt.stats_nhead, 
    #                verbose=opt.verbose)
    # nt occurred per position
    # row: base position
    # col: alphabet order
    LOGGER.timeit(label='preprocess')
    contents, quality_scores_mat, length_dist =  stats_preprocess(
                        seqslist, max_length, alphabet=nucleic_alphabet)
    LOGGER.report(msg='preprocess - Completed in %.2fs.', label='preprocess')

    # 1. GC % for all
    LOGGER.timeit(label='exp1')
    exp1, d1 = gc_content_all(contents, alphabet=nucleic_alphabet, pth = reportdir)
    LOGGER.info("GC content(%) : {:.3f}".format(d1))
    LOGGER.report(msg='GC contents - Completed in %.2fs.', label='exp1')
    info = [['Files', input_bsnames],
            ['Total', n_seqs],
            ['Analyzed', n_seqs - n_removed],
            ["Removed (<{})".format(trim_min_len), n_removed], 
            ['Nucleic', nucleic_alphabet],
            ['GC content', d1]]
    reporter = StatsAnalysisReport(opt.reportdir, opt.reportfile, 
            infoBox = info, backupCount = 5, title = "Fastq Quality Report")

    #reporter.append(exp1)
    
    # 2. Per base sequence quality
    LOGGER.timeit(label='exp2')
    exp2 = sequence_quality_per_base(quality_scores_mat,
                    pth = reportdir)
    LOGGER.report(msg='Per base sequence quality - Completed in %.2fs.', label='exp2')
    reporter.append(exp2)
    
    # 3. Per sequence quality scores
    LOGGER.timeit(label='exp3')
    exp3 = quality_scores_per_sequence(quality_scores_mat,pth = reportdir)
    LOGGER.report(msg='Per sequence quality scores - Completed in %.2fs.', label='exp3')
    reporter.append(exp3)
    
    # 4. Per base sequence content
    LOGGER.timeit(label='exp4')
    exp4 = sequence_content_per_base(contents, pth = reportdir)
    LOGGER.report(msg='Per base sequence content - Completed in %.2fs.', label='exp4')
    reporter.append(exp4)
    
    # 5. Per base GC content
    LOGGER.timeit(label='exp5')
    exp5 = gc_content_per_base(contents, pth = reportdir)
    LOGGER.report(msg='Per base GC content - Completed in %.2fs.', label='exp5')
    reporter.append(exp5)
    
    # 6. Per sequence GC content
    LOGGER.timeit(label='exp6')
    exp6 = gc_content_per_sequence(seqslist, pth = reportdir)
    LOGGER.report(msg='Per sequence GC content - Completed in %.2fs.', label='exp6')
    reporter.append(exp6)
    
    # 7. Per base N content
    LOGGER.timeit(label='exp7')
    exp7 = N_content_per_base(contents, pth = reportdir)
    LOGGER.report(msg='Per base N content - Completed in %.2fs.', label='exp7')
    reporter.append(exp7)
    
    # 8. Sequence Length Distribution
    LOGGER.timeit(label='exp8')
    exp8 = seq_len_dist(length_dist, pth = reportdir)
    LOGGER.report(msg='Sequence Length Distribution - Completed in %.2fs.', label='exp8')
    reporter.append(exp8)
    
    # 9. Sequence Duplication and Overrepresented sequences
    #   d10 is a tuple of representing the most occurred duplicate read. 
    #      d10 = ("AAAA", count, percentage)
    LOGGER.timeit(label='exp9')
    exp9 = sequence_duplicates(seqslist, pth = reportdir)
    LOGGER.report(msg='Sequence Duplication and Overrepresented sequences - Completed in %.2fs.', label='exp9')
    reporter.append(exp9)
    
    # 10. Kmer Content
    # d11 - Kmer ('str', count, percentage)
    # d12 - nt in the head ('str', count, percentage)
    LOGGER.timeit(label='exp10')
    exp10 = kmer_nheads_reads(seqslist, opt.stats_kmer, opt.stats_nhead, 
                                alphabet=None, verbose=False, pth = reportdir)
    LOGGER.report(msg='Kmer Content and nt in head - Completed in %.2fs.', label='exp10')
    reporter.append(exp10)
    reporter.render()

def save2fastq(outputfile, seqslist):
    LOGGER.timeit(label='write')
    with open(outputfile, 'w') as fhandle:
        for seq in seqslist:
            fhandle.write(seq.toFastq())
    LOGGER.report(msg='Write to file - Completed in %.2fs.', label='write')

def parse_cmd(argv):
    """
        Parse command line arguments
    """
    usage = 'usage: %prog [optioins] --fastq xx.fastq' 
    parser = OptionParser(formatter=IndentedHelpFormatterWithNL(), 
                        add_help_option=True, usage=usage, version=__version__)
    # --- 1. general options ---
    general_opts = OptionGroup(parser, "General options")
    general_opts.add_option(
        "--actions", dest="actions", action="append", 
        default=[], 
        help="""Action options define a serials of operations in order specified by users. Actions : 'stats', 'tbs_head', 'tbs_tail', 'tbl_head','tbl_tail'. default is 'stats' only.\n\n  'stats' - generate statistics report in html.\n\n 'tbs_head' - trimming sequences by quality score from heads.\n\n 'tbs_tail' - trimming sequences by quality score from tails.\n\n 'tbl_head' - trimming sequences by length from heads.\n\n 'tbl_tail' - trimming sequences by length from tails.\n\n 'save'     - save current result into a file with fastq format."""
    )
    
    general_opts.add_option(
        "--log", dest='logfilename', default='biTK',
        type='string',
        help="The name of a log file. If not specified, the name will" + 
             " be biTK.log"
    )
    general_opts.add_option(
        "--log-level", dest='log_level', default='info', type='choice', 
        choices=['debug', 'info', 'warnings','error', 'critical','none'], 
        help="The name of a log file. If not specified, the name will" + 
             " be biTK.log"
    )
    general_opts.add_option(
        "-v", "--verbose",
        action="store_true", dest="verbose", default=False,
        help="Show verbose info"
    )
    # --- 1. input and output ---
    inout_opts = OptionGroup(parser, "Input/output settings")
    inout_opts.add_option(
        "--fastqfile", dest="fastqfile", action="append", default=[], 
        help="The input file in fastq format[REQUIRED]."
    )
    inout_opts.add_option(
        "--fastq-type", dest="fastq_type", type='choice',
        choices=['paired', 'single'], default='single',
        help="The type of fastq files: single or paired, default is single."
    )
    inout_opts.add_option(
        "--quality-score-fmt", dest="quality_score_fmt", 
        type='choice', choices=['phred64', 'phred33'], default='phred33',
        help="The format of quality scores used in fastq format. 'phred64' or 'phred33'"
    )
    inout_opts.add_option(
        "-o", "--out", dest="outfile", type='string',
        default='biTK.trimed.fastq',
        help="The output file after trimming sequences."
    )
    # --- 1.1 report setting ---
    report_opts = OptionGroup(parser, "Statistical report settings")
    report_opts.add_option(
        "--reportdir", dest="reportdir", type='string',
        default='./report',
        help="The output path of html report file."
    )
    report_opts.add_option(
        "--reportfile", dest="reportfile", type='string',
        default='biTK_stats_report.html',
        help="The file name of the report."
    )

    # --- 2. trim options ---
    trim_opts = OptionGroup(parser, "Trimming options")
    trim_opts.add_option(
        "--trim-minlen", dest='trim_min_len',
        type='int',
        help='Remove sequences that are shorter than minimum length'
    )
    trim_opts.add_option(
        "--trim-inplace", dest='trim_inplace',
        action='store_true', default=False,
        help='Triming sequences inplace without making a duplicate data'
    )
    # --- 2.1 trim by length ---
    trim_opts.add_option(
        "--trim-len-head", dest='trim_len_head',
        type='int', default=10,
        help='Trim n bases in the head of a sequence, default is 10.'
    )
    trim_opts.add_option(
        "--trim-len-tail", dest='trim_len_tail',
        type='int', default=10,
        help='Trim n bases in the tail of a sequence, default is 10.'
    )
    # --- 2.2 trim by quality score ---
    trim_opts.add_option(
        "--trim-qs-cutoff", dest='trim_quality_score_cutoff',
        type='int', default=30,
        help='Trim n bases in the head of a sequence, default is 30.'
    )
    trim_opts.add_option(
        "--trim-quality-head", dest='trim_quality_head',
        type='int',
        help='Trim bases with quality score of less than threshold in the head of a sequence'
    )
    trim_opts.add_option(
        "--trim-quality-tail", dest='trim_quality_tail',
        type='int',
        help='Trim bases with quality score of less than threshold in the tail of a sequence.'
    )
    trim_opts.add_option(
        "-a", "--trim_adapter", dest='trim_adapter',
        type='string',
        help='Trim an adapter sequence.'
    )

    # --- 3. statistical options ---
    statistics_opts = OptionGroup(
        parser, 
        "Statistics options",
        "Option arguments for statistics."
    )
    statistics_opts.add_option(
        "-k", "--kmer", dest='stats_kmer',
        type='int', default=3,
        help='Occurrencies of k-mer continuous bases'
    )
    statistics_opts.add_option(
        "-n", "--nhead", dest='stats_nhead',
        type='int', default=3,
        help='Occurrencies of k-mer continuous bases'
    )
    # --- 4. Multithreads options ---
    multihtreads_opts = OptionGroup(
        parser, 
        "Multithreads options for loading data",
        "Option arguments for multithreads."
    )
    multihtreads_opts.add_option(
        "--mcpu", dest="mcpu", type='int', default=None,
        help="The number of threads, default is the number of cores."
    )
    multihtreads_opts.add_option(
        "--chks", dest="chks", type='int', default=0,
        help="The chunk size for parallel works, default is 0. The size will be determined automatically. "
    )

    parser.add_option_group(general_opts)
    parser.add_option_group(inout_opts)
    parser.add_option_group(report_opts)
    parser.add_option_group(trim_opts)
    parser.add_option_group(statistics_opts)
    parser.add_option_group(multihtreads_opts)

    options, arguments = parser.parse_args(argv)

    if arguments == 0:
        print ("Error: no arguments found")
        parser.print_help()
        sys.exit(1)
    
    if options.fastq_type not in ['single', 'paired']:
        print ("Error: unknown specified the type of input fastq file. 'single' or 'paired'")
        parser.print_help()
        sys.exit(1)

    # check input fastq files
    if not options.fastqfile:
        print ("Error: do not specify an input fastq file")
        parser.print_help()
        sys.exit(1)
    else:
        #check number
        if (options.fastq_type == 'single') and (len(options.fastqfile)!=1):
            print ("Error: --fastq-type is single but two fastq files were specified.")
            parser.print_help()
            sys.exit(1)
        elif (options.fastq_type == 'paired') and (len(options.fastqfile)!=2):
            print ("Error: --fastq-type is paired but one fastq files was specified.")
            parser.print_help()
            sys.exit(1)
        # check input file existance and readable
        file_test = True
        for f in options.fastqfile:
            if not(isfile(f) and isReadable(f)):
                print("{} does not exist or is unreadable.".format(f))
                file_test = False
        if not file_test:
            sys.exit(1)

    # check actions
    if not options.actions:
        options.actions = ['stats']
    else:
        for act in options.actions:
            if act not in ACTIONS:
                print("Options: unknown specified action. see actions in help")
                parser.print_help()
                sys.exit(1)

    return options

def main(argv):

    global LOGGER
    global ThPool
    # parse command line arguments
    cmds = ' '.join(argv)
    opt = parse_cmd(argv)
    # general opts
    verbose = opt.verbose
    actions = opt.actions
    # input and output
    inputfile = opt.fastqfile[0]
    paired_file = None if opt.fastq_type == 'single' else opt.fastq_type[1]
    outfile = opt.outfile
    fastq_type = opt.fastq_type
    quality_score_fmt = opt.quality_score_fmt
    # input file path, name, and ext
    inputfile_path = os.path.dirname(os.path.abspath(inputfile))
    inputfilename = os.path.basename(inputfile)
    # supported compressed format
    input_ext = os.path.splitext(inputfilename)[-1].lower()[1:]
    compressed_fmt = None
    if input_ext in ['gz', 'bz2']:
        compressed_fmt = input_ext
    inplace = opt.trim_inplace
    # trim opts
    trim_min_len       =  opt.trim_min_len
    trim_quality_head  =  opt.trim_quality_head
    trim_quality_tail  =  opt.trim_quality_tail
    trim_len_head      =  opt.trim_len_head
    # ThreadPool
    mcpu = opt.mcpu if opt.mcpu else cpu_count()
    chks = opt.chks
    # io_plug
    if mcpu == 1:
        io_plug = 'FastQIO'
    else:
        io_plug = 'FastQIO_multithread'

    # Logger settings
    logfile = opt.logfilename
    log_label = inputfilename
    LOGGER = Console(log_label, prefix='@>', console=opt.log_level)
    #LOGGER = Console(log_label, prefix='@>')
    LOGGER.start(logfile)

    LOGGER.info(' '.join(argv))
    LOGGER.timeit(label='Load fastq')
    seqslist, removed_list, max_length = load_data(inputfile, trim_min_len,
                        io_plug = io_plug, 
                        fastq_type = fastq_type,
                        alphabet=nucleic_alphabet,
                        quality_score_fmt=quality_score_fmt,
                        compressed=compressed_fmt, 
                        nthreads = mcpu,
                        chks = chks, 
                        verbose=verbose)
    if paired_file:
        seqslist_p, removed_list_p, max_length_p = load_data(inputfile, trim_min_len,
                        io_plug = io_plug, 
                        fastq_type = fastq_type,
                        alphabet=nucleic_alphabet,
                        quality_score_fmt=quality_score_fmt,
                        compressed=compressed_fmt, 
                        nthreads = mcpu,
                        chks = chks, 
                        verbose=verbose)
    LOGGER.report(msg='Load fastq Completed in %.2fs.', label='Load fastq')

    len_actions = len(opt.actions)
    seqslist_new = None
    if len_actions == 1 and opt.actions[0] == 'stats':
        statistics_report(seqslist, removed_list, max_length, cmds, opt)
    else:
        for act in opt.actions:
            if act == 'stats':
                statistics_report(seqslist, removed_list, max_length, cmds, opt)
            elif act == 'tbs_head':
                if inplace:
                    tbs_head(seqslist, trim_quality_head, 
                                inPlace=inplace, verbose=verbose)
                else:
                    seqslist_new = tbs_head(seqslist, trim_quality_head, 
                                inPlace=inplace, verbose=verbose)
            elif act == 'tbs_tail':
                if inplace:
                    tbs_tail(seqslist, trim_quality_tail, 
                                inPlace=inplace, verbose=verbose)
                else:
                    seqslist_new = tbs_tail(seqslist, trim_quality_tail, 
                                inPlace=inplace, verbose=verbose)
            elif act == 'tbl_head':
                if inplace:
                    tbl_head(seqslist, trim_len_head, 
                                inPlace=inplace, verbose=verbose)
                else:
                    seqslist_new = tbl_head(seqslist, trim_len_head, 
                                inPlace=inplace, verbose=verbose)
            elif act == 'tbl_tail':
                if inplace:
                    tbl_tail(seqslist, trim_len_tail, 
                            inPlace=inplace, verbose=verbose)
                else:
                    seqslist_new = tbl_tail(seqslist, trim_len_tail, 
                            inPlace=inplace, verbose=verbose)
            elif act == 'save':
                if inplace:
                    save2fastq(outfile, seqslist)
                else:
                    save2fastq(outfile, seqslist_new)
    # Shutdown ThreadPool
#    ThPool.shutdown()

if __name__ == '__main__':
    main(sys.argv)
