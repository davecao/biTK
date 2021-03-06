# -*- coding: utf-8 -*-
from __future__ import print_function

from functools import wraps, partial
import os
import sys
import re
import pprint
import time
import mmap
from itertools import islice
import biTK
from biTK.ngs.sequence import Sequence, Alphabet, nucleic_alphabet,\
                            PHRED_ALPHABET
from biTK import PY3K
from biTK.ngs.io.pathtools import OPEN
from biTK.ngs.utils import grouper
# from biTK.ngs.utils import profile

from multiprocessing import cpu_count
# from multiprocessing.pool import ThreadPool
# memory
# from memory_profiler import profile

# joblib
# import joblib.parallel
from collections import defaultdict
from tempfile import mkdtemp
from joblib import Parallel, delayed, Memory

# concurrent futures ProcessPoolExecutor
import concurrent.futures
from signal import signal, SIGPIPE, SIG_DFL
# reset
signal(SIGPIPE, SIG_DFL)

if PY3K:
    # On Python 3, this will be a unicode StringIO
    import io
    from io import StringIO
    StringTypes = bytes
    FileType = io.IOBase
    ClassType = type
    InstanceType = object
else:
    # On Python 2 this will be a (bytes) string based handle.
    # Note this doesn't work as it is unicode based:
    # from io import StringIO
    try:
        from cStringIO import StringIO
    except ImportError:
        from StringIO import StringIO
    from types import StringTypes, FileType, ClassType, InstanceType

cachedir = mkdtemp()
memory = Memory(cachedir=cachedir, mmap_mode='r', verbose=0)

# __all__ = [
#        "AlignIO", "MultiFastaIO", "ClustalWIO",
#        "FastQIO", FastQIO_multithread"]

__all__ = ["AlignIO", "sequenceBuilder"]


def stdrepr(obj,  attributes=None, name=None):
    """Create a standard representation of an object."""
    if name is None:
        name = obj.__class__.__name__
    if attributes is None:
        attributes = obj.__class__.__slots__
    args = []
    for a in attributes:
        if a[0] == '_':
            continue
        args.append('%s=%s' % (a, repr(getattr(obj, a))))
    args = ',\n'.join(args).replace('\n', '\n    ')
    return '%s(\n    %s\n)' % (name, args)


# decorator borrowed from Mozilla mxr
def abstractmethod(method):
    if PY3K:
        line = method.__code__.co_firstlineno
        filename = method.__code__.co_filename
    else:
        line = method.func_code.co_firstlineno
        filename = method.func_code.co_filename

    @wraps(method)
    def not_implemented(*args, **kwargs):
        raise NotImplementedError(
            'Abstract method %s at File "%s", line %s'
            'should be implemented by a concrete class' %
            (repr(method), filename, line))
    return not_implemented


class ClassRegistry(type):
    """ Register all subclasses """
    def __init__(cls, name, bases, nmspc):
        super(ClassRegistry, cls).__init__(name, bases, nmspc)
        if not hasattr(cls, 'registry'):
            cls.registry = set()
        cls.registry.add(cls)
        cls.registry -= set(bases)  # Remove base classes

    # Meta methods, called on class objects:
    def __iter__(cls):
        return iter(cls.registry)

    def __str__(cls):
        if cls in cls.registry:
            return cls.__name__
        return cls.__name__ + ":" + ', '.join([sc.__name__ for sc in cls])


class AlignIO(object):
    """
    .. note::
        This class uses two patterns, composite and registry

        supported format:
            FASTA - MultiFastaIO,
            clustalW - ClustalWIO,
            Phylip - PhylipIO
            FASTQ - FastQIO

    .. ipython:: python
        from biTK.ngs.sequence import Alphabet, generic_alphabet, protein_alphabet, nucleic_alphabet
        from biTK.ngs.io import AlignIO
        filename="./Desktop/multifasta.aln"
        parser=AlignIO(filename, ConcreteIO='ClustalWIO')
        seqs=parser.parse(alphabet=protein_alphabet, isAligned=True)

    """
    __metaclass__ = ClassRegistry

    def __init__(self, *args, **kwargs):
        """Initialization.

        Args:
           IO (class):

        Kwargs:
        """
        super(AlignIO, self).__init__()
        ConcreteIO_cls = kwargs.pop("ConcreteIO", None)
        if ConcreteIO_cls is None:
            return
#        if isinstance(ConcreteIO, ClassType):
#            # developped by a user or subclasses of Kernel
#            if ConcreteIO in self.registry:
#                self.__ConcreteIO = ConcreteIO(*args, **kwargs)
#        elif isinstance(ConcreteIO, InstanceType):

#            self.__ConcreteIO = ConcreteIO
#        elif self.__isstr(ConcreteIO):
#            # is string
        self.__ConcreteIO = self.__create(ConcreteIO_cls, *args, **kwargs)

    def __create(self, clsname, *args, **kwargs):
        """ Create an instance for a given class in str or name
        """
        obj = None
        for cls in self.registry:
            if clsname == cls.__name__ or clsname == cls:
                obj = cls(*args, **kwargs)
        if obj:
            # Initialize object
            obj.__init__(*args, **kwargs)
        else:
            print("Unknown class {}".format(clsname))
        return obj

    def __str__(self):
        IO_func_str_ = self.__getattr__(self.__ConcreteIO, '__str__')
        return IO_func_str_()

    def __repr__(self):
        IO_func_str_ = self.__getattr__(self.__ConcreteIO, '__repr__')
        return IO_func_str_()

    def __call__(self, *args, **kwargs):
        return self.__ConcreteIO(*args, **kwargs)

    def __getattr__(self, attr):
        """ get delegation to the object """
        try:
            return self.__ConcreteIO.__getattribute__(attr)
        except AttributeError:
            raise AttributeError(
                '{0} object has no attribute `{1}`'
                .format(self.__class__.__name__, attr))

    def __isstr(self, s):
        try:
            return isinstance(s, basestring)
        except NameError:
            return isinstance(s, str)


class IOBase(object):
    """ Base class for Conrete IO classes """
    def __init__(self, handle, *args, **kwargs):
        """Initialization.
        Args:
           *args:
        Kwargs:
            compressed : zip | bz2 | gzip
        """
        super(IOBase, self).__init__()
        compressed = kwargs.get('compressed', None)
        if compressed:
            fhandle = StringIO()
            with OPEN[compressed](handle, 'r') as f:
                fhandle.write(f.read())
            handle = fhandle

        if isinstance(handle, FileType):
            # file handle
            self.handle = handle
        elif isinstance(handle, StringTypes):
            # is string
            try:
                self.handle = open(handle, 'r')
            except IOError as e:
                errno, strerr = e.args
                print ("I/O error({0}): {1}".format(errno, strerr))
            except:
                print("Unexpected error:{0}".format(sys.exc_info()[0]))
        elif hasattr(handle, "read"):
            self.handle = handle
        else:
            raise IOError("Unknown argument: a file handle or a string")
        return self

    @abstractmethod
    def parse(self, alphabet=None, isAligned=False):
        """
        Args:
            handle   -- handle to the file, or the filename as a string
            alphabet -- The expected alphabet of the data, if given

        Returns:
            seq -- a list of bilab.sequence.Sequence
        """


class MultiFastaIO(IOBase, AlignIO):
    """ Derived class
    To read multiple alignments in Fasta format
    e.g.
        >id1:...
        ACDEF....
        DDD--...
        >id2:...
        -CD-F...
        DD-D-...
    """
    def __init__(self, handle):
        super(MultiFastaIO, self).__init__(handle)

    def parse(self, alphabet=None, isAligned=False):
        """ Impletementation of the abstract method defined in IOBase
        Args:
            handle   -- handle to the file, or the filename as a string
            alphabet -- The expected alphabet of the data, if given

        Returns:
            seq -- a list of bilab.sequence.Sequence
        """
        handle = self.handle
        alphabet = Alphabet(alphabet)

        def build_seq(seq, alphabet, header, header_lineno,
                      comments, isAligned=False):
            try:
                # Create a bilab.sequence object
                title_fields = re.split("\||:", header)
                seq_id = title_fields[0]
                if comments:
                    header += '\n' + '\n'.join(comments)
                s = Sequence(seq, alphabet, name=seq_id, description=header)
            except ValueError:
                raise ValueError(
                    "Failed to parse the file at the line %d: "
                    "Character not in alphabet:%s" % (header_lineno, alphabet))
            return s

        # loop file handle
        if not isinstance(handle, file):
            print("IOError: argument is not a file handle")
            sys.exit(0)
        seqs = []
        seq_str = ""
        comments = []
        header = None
        header_lineno = 0
        for lineno, line in enumerate(handle):
            line = line.strip("\n")
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    # Create Sequence object
                    s = build_seq(seq_str, alphabet, header,
                                  header_lineno, comments, isAligned=isAligned)
                    seqs.append(s)
                    seq_str = ""
                    header = None
                    comments = []
                header = line[1:]
                header_lineno = lineno
            elif line.startswith(';'):
                comments.append(line[1:])
            else:
                seq_str += line
        # Store last one
        s = build_seq(seq_str, alphabet, header,
                      header_lineno, comments, isAligned=isAligned)
        seqs.append(s)
        return seqs


class ClustalWIO(IOBase, AlignIO):
    """
    Derived class
    To read multiple alignments in clustalw format
    e.g.
CLUSTAL W (1.81) multiple sequence alignment


CXCR3_MOUSE       --------------------------LENSTSPYDYGENESD-------FSDSPPCPQDF
BLR_HUMAN         --------------------------LENLEDLF-WELDRLD------NYNDTSLVENH-
CXCR1_HUMAN       --------------------------MSNITDPQMWDFDDLN-------FTGMPPADEDY
CXCR4_MURINE      -----------------------------------YTSDN---------YSGSGDYDSNK
                                                     :  :          :..     ..

CXCR3_MOUSE       -SL-------NFDRTFLPALYSLLFLLGLLGNGAVAAVLLSQRTALSSTDTFLLHLAVAD
BLR_HUMAN         --LC-PATMASFKAVFVPVAYSLIFLLGVIGNVLVLVILERHRQTRSSTETFLFHLAVAD
CXCR1_HUMAN       -SPC-MLETETLNKYVVIIAYALVFLLSLLGNSLVMLVILYSRVGRSVTDVYLLNLALAD
CXCR4_MURINE      -EPC-RDENVHFNRIFLPTIYFIIFLTGIVGNGLVILVMGYQKKLRSMTDKYRLHLSVAD
                             :.  .:   * ::** .::**  *  ::   :   * *: : ::*::**

CXCR3_MOUSE       VLLVLTLPLWAVDAA-VQWVFGPGLCKVAGALFNINFYAGAFLLACISFDRYLSIVHATQ
BLR_HUMAN         LLLVFILPFAVAEGS-VGWVLGTFLCKTVIALHKVNFYCSSLLLACIAVDRYLAIVHAVH
CXCR1_HUMAN       LLFALTLPIWAASKV-NGWIFGTFLCKVVSLLKEVNFYSGILLLACISVDRYLAIVHATR
CXCR4_MURINE      LLFVITLPFWAVDAM-ADWYFGKFLCKAVHIIYTVNLYSSVLILAFISLDRYLAIVHATN
                  :*:.: **: ...     * :*  ***..  :  :*:*.. ::** *:.****:****..

    """
    class Token(object):
        """Represents the items returned by a file scanner, normally processed
        by a parser.

        Attributes :
        o typeof    -- a string describing the kind of token
        o data      -- the value of the token
        o lineno    -- the line of the file on which the data was found
                       (if known)
        o offset    -- the offset of the data within the line (if known)
        """
        __slots__ = ['typeof', 'data', 'lineno', 'offset']

        def __init__(self, typeof, data=None, lineno=-1, offset=-1):
            self.typeof = typeof
            self.data = data
            self.lineno = lineno
            self.offset = offset

        def __repr__(self):
            return stdrepr(self)

        def __str__(self):
            coord = str(self.lineno)
            if self.offset != -1:
                coord += ':'+str(self.offset)
            coord = coord.ljust(7)
            return (coord + '  ' + self.typeof + ' : ').ljust(32) +\
                    str(self.data or '')

    def __init__(self, handle):
        super(ClustalWIO, self).__init__(handle)

    def parse(self, alphabet=None, isAligned=False):
        """ Impletementation of the abstract method defined in IOBase
            referred to clustal_io.py in weblogo
        Args:
            handle   -- handle to the file, or the filename as a string
            alphabet -- The expected alphabet of the data, if given

        Returns:
            seq -- a list of bilab.sequence.Sequence
        """

        def scan(handle):
            """Scan a clustal format MSA file and yield tokens.
            The basic file structure is

            begin_document
                header?
               (begin_block
                   (seq_id seq seq_index?)+
                   match_line?
               end_block)*
            end_document

            Usage:
                for token in scan(clustal_file):
                    do_something(token)
            """
            header, body, block = range(3)
            yield self.Token("begin")
            leader_width = -1
            state = header
            for L, line in enumerate(handle):
                if state == header:
                    if line.isspace():
                        continue
                    m = header_line.match(line)
                    state = body
                    if m is not None:
                        yield self.Token("header", m.group())
                        continue
                    # Just keep going and hope for the best.
                    # else :
                        # raise ValueError("Cannot find required header")

                if state == body:
                    if line.isspace():
                        continue
                    yield self.Token("begin_block")
                    state = block
                    # fall through to block

                if state == block:
                    if line.isspace():
                        yield self.Token("end_block")
                        state = body
                        continue

                    m = match_line.match(line)
                    if m is not None:
                        yield self.Token("match_line", line[leader_width:-1])
                        continue

                    m = seq_line.match(line)
                    if m is None:
                        raise ValueError(
                            "Parse error on line: %d (%s)" % (L, line))
                    leader_width = len(m.group(1))
                    yield self.Token("seq_id", m.group(1).strip())
                    yield self.Token("seq", m.group(2).strip())
                    if m.group(3):
                        yield self.Token("seq_num", m.group(3))
                    continue

                # END state blocks.
                # If I ever get here something has gone terrible wrong
                raise RuntimeError()

            if state == block:
                yield self.Token("end_block")
            yield self.Token("end")
            return

        handle = self.handle
        header_line = re.compile(r'(CLUSTAL.*)$')
        # (sequence_id) (Sequence) (Optional sequence number)
        seq_line = re.compile(r'(\s*\S+\s+)(\S+)\s*(\d*)\s*$')
        # Saved group includes variable length leading space.
        # Must consult a seq_line to figure out how long the leading space
        # is since the maximum CLUSTAL ids length (normally 10 characters)
        # can be changed.
        match_line = re.compile(r'([\s:\.\*]*)$')
        alphabet = Alphabet(alphabet)
        seq_ids = []
        seqs = []
        block_counts = 0
        data_len = 0
        for token in scan(handle):
            if token.typeof == "begin_block":
                block_count = 0
            elif token.typeof == "seq_id":
                if len(seqs) <= block_count:
                    seq_ids.append(token.data)
                    seqs.append([])
            elif token.typeof == "seq":
                if not alphabet.alphabetic(token.data):
                    raise ValueError(
                        "Character on line: %d not in alphabet: %s : %s" % (
                        token.lineno, alphabet, token.data))
                seqs[block_count].append(token.data)
                if block_count == 0:
                    data_len = len(token.data)
                elif data_len != len(token.data):
                    raise ValueError("Inconsistent line lengths")

                block_counts += 1
        seqs = [
            Sequence("".join(s), alphabet, name=i)
            for s, i in zip(seqs, seq_ids)]
        return seqs


class PhylipIO(IOBase, AlignIO):
    """ Derived class
    To read multiple alignments in PHYLIP format
1. First line contains number of species and number of characters in a species
   sequence.
2. A user tree may appear at end of the file
3. options
  U - User Tree
  G - Global
  J - Jumble
  O - Outgroup
  T - Threshold
  M - Multiple data sets
  W - Weights
-------------------------------------------
  6   50   W
W         0101001111 0101110101 01011
dmras1    GTCGTCGTTG GACCTGGAGG CGTGG
hschras   GTGGTGGTGG GCGCCGGCCG TGTGG
ddrasa    GTTATTGTTG GTGGTGGTGG TGTCG
spras     GTAGTTGTAG GAGATGGTGG TGTTG
scras1    GTAGTTGTCG GTGGAGGTGG CGTTG
scras2    GTCGTCGTTG GTGGTGGTGG TGTTG

0101001111 0101110101 01011
GTCGTCGTTG GACCTGGAGG CGTGG
GTGGTGGTGG GCGCCGGCCG TGTGG
GTTATTGTTG GTGGTGGTGG TGTCG
GTAGTTGTAG GAGATGGTGG TGTTG
GTAGTTGTCG GTGGAGGTGG CGTTG
GTCGTCGTTG GTGGTGGTGG TGTTG

1
((dmras1,ddrasa),((hschras,spras),(scras1,scras2)));
-----------------------------------------------------
    """
    def __init__(self, handle):
        super(PhylipIO, self).__init__(handle)

    def parse(self, alphabet=None, isAligned=False):
        """ Parse Phylip format, return a list of bilab.sequence.Sequence
            See weblog's phylip_io
        Args:

        Kwargs:
            alphabet - bilab.sequence.Alphabet

        return:
            a list of bilab.sequence.Sequence
        """
        seqs = []
        idents = []
        num_seq = 0
        num_total_seq = 0  # length of sequence of 1 species
        tracker = 0  # track what sequence the line is on
        usertree_tracker = 0  # track usertree lines
        options = ''  # options
        num_options = 0  # number/lens of options - U
        handle = self.handle
        for lineno, currline in enumerate(handle):
            # split into fields
            line = currline.strip("\n").split()
            if line == []:
                continue

            if (line[0].isdigit() and len(line) == 1 and
                    len(seqs) == num_seq and
                    len(seqs[0]) == num_total_seq):
                usertree_tracker = int(line[0])
                pass  # identifies usertree
            elif num_options > 0:
                if len(seqs) < num_seq:
                    if line[0][0] in options:
                        num_options -= 1
                        pass
                    else:
                        raise ValueError("Not an option, but it should be one")
                else:
                    num_options -= 1
                    pass
            elif usertree_tracker > 0:
                # skip usertree basically
                if len(seqs[num_seq - 1]) == num_total_seq:
                    usertree_tracker -= 1
                    pass
                else:
                    raise ValueError('User Tree in Wrong Place')
            elif line[0].isdigit():
                if len(line) >= 2 and len(seqs) == 0:
                    # identifies first line
                    # number of sequences
                    num_seq = int(line[0])
                    # length of sequences
                    num_total_seq = int(line[1])
                    if len(line) > 2:
                        options = (''.join(line[2:]))
                        num_options = len(options) - options.count('U')
                    # else:
                    #    raise ValueError('parse error')
            elif num_options == 0:
                if num_seq == 0:
                    raise ValueError("Empty File, or possibly wrong file")
                elif tracker < num_seq:
                    if num_seq > len(seqs):
                        seqs.append(''.join(currline[10:].split()))
                        # removes species name
                        idents.append(currline[0:10].strip())
                        tracker += 1
                    else:
                        seqs[tracker] += (''.join(line))
                        tracker += 1
                    if tracker == num_seq:
                        tracker = 0
                        num_options = len(options) - options.count('U')
        if len(seqs) != len(idents) or len(seqs) != num_seq:
            raise ValueError("Number of different sequences wrong")
        sequence = []
        for i in range(0, len(idents)):
            if len(seqs[i]) == num_total_seq:
                sequence.append(
                    Sequence(seqs[i], alphabet=alphabet,
                             name=idents[i], isAligned=isAligned))
            else:
                raise ValueError("extra sequence in list")
        return sequence


class FastQIO(IOBase, AlignIO):
    """ Derived class
    To read sequences in FastQ format
    e.g.
        @id1...
        ACGT....
        +id1....
        !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
        @id2:...
        CAGT...
        +id2...
        !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
        ...
        ...
    """
    def __init__(self, handle, *args, **kwargs):
        super(FastQIO, self).__init__(handle, *args, **kwargs)

#    @profile
    def parse(self, alphabet=nucleic_alphabet,
              quality_score_fmt='phred33', **kwargs):
        """ Impletementation of the abstract method defined in IOBase
        Args:
            handle   -- handle to the file, or the filename as a string
            alphabet -- The expected alphabet of the data, if given

        Returns:
            seq -- a list of bilab.ngs.Sequence
        """
        # from itertools import islice
        handle = self.handle
        alphabet = Alphabet(alphabet)
        q_alphabet = PHRED_ALPHABET[quality_score_fmt]
        sequence_nlimit = kwargs.get('sequence_nlimit', 200000)

#        def build_seq(raw_seq, quality_seq, header, option_id):
#            try:
#                # Create a bilab.ngs.sequence object
#                title_fields = re.split("[\s|:]", header)
#                sequence_id = title_fields[0]
#
#                s = Sequence(raw_seq, quality_seq,
#                                alphabet=alphabet,
#                                quality_alphabet = q_alphabet,
#                                name=sequence_id,
#                                optionID=option_id,
#                                quality_format=quality_score_fmt,
#                                description=header)
#            except ValueError:
#                raise ValueError(
#                    "Failed to parse the file at the line %d: "
#                    "Character not in alphabet: %s"
#                    " %s"%(lineno, alphabet, raw_seq))
#            return s

        # loop file handle
        if not isinstance(handle, file):
            if not hasattr(self.handle, "read"):
                print("IOError: argument is not a recognized handle")
                sys.exit(0)
            else:
                # has read method -- StringIO
                handle = self.handle.getvalue().split('\n')
#        seqs = []
#        raw_sequence = ""
#        quality_sequence = ""
#        option_id = ""
#        state = -1
#        header = None
#        header_lineno = 0
#        id_, oid_, quality = range(3)
#        uniq_name = re.split('[:|.\s]', handle.readline())[0]
#        handle.seek(0)
#        if sequence_nlimit <0:
#            handle_lines = handle
#        else:
#            handle_lines = islice(handle, 0, sequence_nlimit*4)
#        for lineno, line in enumerate(handle_lines):
#            line = line.strip("\n")
#            if not line:
#                continue
#            if line.startswith(uniq_name):
#                state = id_
#                if header is not None:
#                    # Create Sequence object
#                    s = build_seq(raw_sequence, quality_sequence, header,
#                                    option_id)
#                    seqs.append(s)
#                    raw_sequence = ""
#                    quality_sequence = ""
#                    header =None
#                header = line[1:]
#                header_lineno = lineno
#            elif line.startswith('+'):
#                option_id = line[1:]
#                state = quality
#            #elif line.startswith('!'):
#            #    quality_sequence = line[1:]
#            #    state = quality
#            else:
#                if state == id_:
#                    raw_sequence += line
#                elif state == quality:
#                    quality_sequence += line
#
#        # Store last one
#        s = build_seq(raw_sequence, quality_sequence,header,
#                                    option_id)
#        seqs.append(s)
        seqs = []
        seqs_append = seqs.append
        # m_handle = mmap.mmap()
        if sequence_nlimit <0:
            handle_lines = handle
        else:
            handle_lines = islice(handle, 0, sequence_nlimit*4)

        for header, raw_seq, option_id, quality_seq in grouper(4, handle_lines):
            try:
                if header[0] != "@" or option_id[0] != "+":
                    print("Malformated fastq file:")
                    print("{}\n{}\n{}\n\{}".format(
                        header, raw_seq, option_id, quality_seq))
                    sys.exit(1)

                title_fields = re.split("[\s|:]", header)
                sequence_id = title_fields[0]
                s_obj = Sequence(
                            raw_seq.strip(),
                            quality_seq.strip(),
                            alphabet=alphabet,
                            name=sequence_id,
                            optionID=option_id,
                            qs_fmt=quality_score_fmt,
                            description=header.strip())
                seqs_append(s_obj)
            except ValueError:
                raise ValueError(
                    "Character not in alphabet: %s %s" % (alphabet, raw_seq))
        handle.close()
        return seqs


def sequenceBuilder(header, raw_seq, option_id, quality_seq, alphabet,
                    quality_score_fmt):
    try:
        # Create a bilab.ngs.sequence object
        title_fields = re.split("[\s|:]", header)
        sequence_id = title_fields[0]

        # alphabet = kwargs.get('alphabet', nucleic_alphabet)
        # quality_score_fmt = kwargs.get('quality_score_fmt', 'phred33')
        s_obj = Sequence(
                    raw_seq.strip(),
                    quality_seq.strip(),
                    alphabet=alphabet,
                    name=sequence_id,
                    optionID=option_id,
                    qs_fmt=quality_score_fmt,
                    description=header.strip())
    except ValueError:
       raise ValueError(
            "Character not in alphabet: %s %s" % (alphabet, raw_seq))
    return s_obj


class JoblibCallBack(object):
    completed = defaultdict(int)

    def __init__(self, index, parallel):
        self.index = index
        self.parallel = parallel

    def __call__(self, index):
        JoblibCallBack.completed[self.parallel] += 1
        print("done with {}".format(JoblibCallBack.completed[self.parallel]))
        if self.parallel._original_iterable:
            self.parallel.dispatch_next()


class FastQIO_multithread(IOBase, AlignIO):
    """ Derived class
    To read sequences in FastQ format
    e.g.
        @id1...
        ACGT....
        +id1....
        !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
        @id2:...
        CAGT...
        +id2...
        !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
        ...
        ...
    """
    def __init__(self, handle, *args, **kwargs):
        super(FastQIO_multithread, self).__init__(handle, *args, **kwargs)

    def parse(self, alphabet=nucleic_alphabet, quality_score_fmt='phred64',
              sequence_nlimit=200000,
              nthreads=None, chunksize=0, verbose=0):
        """ Impletementation of the abstract method defined in IOBase
        Args:
            handle   -- handle to the file, or the filename as a string
            alphabet -- The expected alphabet of the data, if given

        Returns:
            seq -- a list of bilab.ngs.Sequence
        """
        handle = self.handle
        alphabet = Alphabet(alphabet)
        # ThreadPool: use the number of cores available
        if nthreads is None:
            nthreads = cpu_count()*2
#            pool = ThreadPool(cpu_count()*2)
        else:
            nthreads = nthreads
#            pool = ThreadPool(nthreads)
        if chunksize is None:
            chunksize = 'auto'
        elif isinstance(chunksize, int):
            if chunksize <= 0:
                chunksize = 'auto'

        # loop file handle
        if not isinstance(handle, file):
            if not hasattr(self.handle, "read"):
                print("IOError: argument is not a recognized handle")
                sys.exit(0)
            else:
                # has read method -- StringIO
                handle = self.handle.getvalue().split('\n')
        if sequence_nlimit < 0:
            handle_lines = handle
        else:
            handle_lines = islice(handle, 0, sequence_nlimit*4)
        # using cache makes computations more slower
        # sequenceBuilderCached = memory.cache(sequenceBuilder)
        # joblib.parallel.CallBack = JoblibCallBack
        p = Parallel(n_jobs=nthreads, backend="multiprocessing",
                     batch_size=chunksize, verbose=verbose,
                     temp_folder="/tmp/biTK",
                     max_nbytes='100M', mmap_mode='r')
        func = delayed(sequenceBuilder, check_pickle=True)
        seqs = p(func(
                    header, raw_seq, option_id, quality_seq,
                    alphabet, quality_score_fmt
                    )
                    for header, raw_seq, option_id, quality_seq,
                        alphabet, quality_score_fmt in grouper(
                                4, handle_lines,
                                opts=(alphabet, quality_score_fmt)
                        ))
        handle.close()
        return seqs
