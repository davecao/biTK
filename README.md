biTK package
=============

Prerequisites
-------------

1. Jinja2
2. Matplotlib

Download and Installation
-------------------------
  pysetup run build install_dist --prefix=/path/to/install
or
  python setup.py install


Command Line Options
--------------------

python -m 'biTK.ngs.apps.statsApp'

Usage: statsApp.py [optioins] --fastq xx.fastq

Options: 
  --version             show program's version number and exit 
  -h, --help            show this help message and exit 

  General options: 
    --action=ACTION     The input file in fastq format[REQUIRED]. 
    --fastqfile=FASTQFILE 
                        The input file in fastq format[REQUIRED]. 
    -o OUTFILE, --out=OUTFILE 
                        The output file after trimming sequences. 
    --reportfile=REPORTFILE 
                        The output path of html report file. 
    --quality-score-fmt=QUALITY_SCORE_FMT
                        The format of quality scores used in fastq format.
    --log=LOGFILENAME   The name of a log file. If not specified, the name
                        will be biTK.log
    -v, --verbose       Show verbose info

  Trimming options:
    --trim-minlen=TRIM_MIN_LEN
                        Remove sequences that are shorter than minimum length
    --trim-inplace      Triming sequences inplace without making a duplicate
                        data
    --trim-qs-cutoff=TRIM_QUALITY_SCORE_CUTOFF
                        Trim n bases in the head of a sequence, default is 30.
    --trim-len-head=TRIM_LEN_HEAD
                        Trim n bases in the head of a sequence, default is 10.
    --trim-len-tail=TRIM_LEN_TAIL
                        Trim n bases in the tail of a sequence, default is 10.
    --trim-quality-head=TRIM_QUALITY_HEAD
                        Trim bases with quality score of less than threshold
                        in the head of a sequence
    --trim-quality-tail=TRIM_QUALITY_TAIL
                        Trim bases with quality score of less than threshold
                        in the tail of a sequence.
    -a TRIM_ADAPTER, --trim_adapter=TRIM_ADAPTER
                        Trim an adapter sequence.

  Statistics options:
    Option arguments for statistics.

    -k STATS_KMER, --kmer=STATS_KMER
                        Occurrencies of k-mer continuous bases
    -n STATS_NHEAD, --nhead=STATS_NHEAD
                        Occurrencies of k-mer continuous bases

