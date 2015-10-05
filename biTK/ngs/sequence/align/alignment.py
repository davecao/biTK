# -*- coding: utf-8 -*-

import datetime
import sys

from biTK.ngs.SubMatrice.MatInfo import SubstitutionMatrix as SUBMat
from biTK.ngs.sequence.align import AlignInfo

__all__ = ['Alignment']

class ClassRegistry(type):
    """ Register all subclasses """
    def __init__(cls, name, bases, nmspc):
        super(ClassRegistry, cls).__init__(name, bases, nmspc)
        if not hasattr(cls, 'registry'):
            cls.registry = set()
        cls.registry.add(cls)
        cls.registry -= set(bases) #Remove base classes
    # Meta methods, called on class objects:
    def __iter__(cls):
        return iter(cls.registry)

    def __str__(cls):
        if cls in cls.registry:
            return cls.__name__
        return cls.__name__ + ":" + ', '.join([sc.__name__ for sc in cls])

class Cell(object):
    """
        Store an element of a matrix
        traceback direction:
        diag = 0
        up   = 1
        left = 2
    """
    __slot__ = ['row','col','query_aa','target_aa','direction','value','parent']

    def __init__(self, row, col, target_aa, query_aa, direction, parent, value):
        super(Cell, self).__init__()
        self.col = col
        self.row = row
        self.query_aa = query_aa
        self.target_aa = target_aa
        self.direction = direction
        self.value = value
        self.parent = parent


class DynamicProgMatrix(object):
    """
        The matrix of Dynamic Programming 
                 A         C        A  (Query)
             ____0_________0________0______
            | dia    1|        2|        3|
Target C  0 |  (i-1,  | (i-1,j) |         |
            |   j-1)  |   up    |         |
            |1________|1________|1________|
            | left   1|        2|        3| 
       C  0 | (i,j-1) | (i,j)   |         | 
            |2________|2________|2________|
            |        1|        2|        3|
       A  0 |         |         |         |
            |3________|3________|3________|
        
        currentCell.col = 1
        currentCell.row = 1
        currentCell.query_aa = C
        currentCell.target_aa = C
        currentCell.value = 
        currentCell.parent = 
    Needleman-Wunsch scheme
      qdiag = C(i−1,j−1) + S(i,j)  diag
      qup   = C(i−1,j)   + g(i-1)  up
      qleft = C(i,j−1)   + g(j-1)  left

    """
    def __init__(self, sub_matrix, target, query, align_type, 
                        gap_open, gap_extension, gap_symbol='-'):
        self.sub_mat = sub_matrix
        self.target = target
        self.query = query
        self.align_type = align_type
        self.gap_open = gap_open
        self.gap_extension = gap_extension
        self.gap_symbol = gap_symbol
        self.best_score = -9999.0
        self.best_score_inx = (0,0)
        self.ncols = len(query) + 1
        self.nrows = len(target) + 1
        self.dpMat = {}
        #defaultdict((,))

    def _argmax(self, slist):
        """ 
            Return max value and its index in a list

        """
        #slist = [score_diag, score_up, score_left]
        max_score = max(slist)
        inx = [i for i, j in enumerate(slist) if j==max_score]
        return max_score, inx[0]
    
    def fill(self, strategy='global'):
        """
            Actions for filling dp matrix
        Args:

        Kwargs:
            strategy (str): global or Local
        """
        if strategy == 'global':
            self.fill_global()
        elif strategy == 'local':
            self.fill_local()

    def fill_global(self):
        ncols = self.ncols #len(self.query)
        nrows = self.nrows #len(self.target)
        diag, up, left = range(3)

        self.dpMat[(0, 0)] = Cell(0,0,None, None, None, None, 0.0)
        for j in range(1, ncols):
            # pad cols
            query_aa = self.query[j-1].upper()
            val = self.affine_gap_penalty(j)
            self.dpMat[(0, j)] = Cell(0,j,"", query_aa, left, (0,j-1), val)
        for i in range(1,nrows):
            # pad rows
            target_aa = self.target[i-1].upper()
            val = self.affine_gap_penalty(i)
            self.dpMat[(i, 0)] = Cell(i,0, target_aa, "", up, (i-1, 0), val)

        #Fill the score matrix
        for i in range(1, nrows):
            target_aa = self.target[i-1].upper()
            for j in range(1, ncols):
                query_aa = self.query[j-1].upper()
                # no gap: diag
                score_diag = self.dpMat[(i-1, j-1)].value + \
                             self.sub_mat.get_value(target_aa, query_aa)
                # up
                score_up = self.dpMat[(i-1, j)].value + \
                            self.affine_gap_penalty(i-1)
                # left
                score_left = self.dpMat[(i, j-1)].value + \
                            self.affine_gap_penalty(j-1)

                max_score, direct = self._argmax([score_diag, 
                                                score_up, score_left])

                if max_score > self.best_score:
                    self.best_score = max_score
                    self.best_score_inx = (i,j)

                if direct == diag:
                    self.dpMat[(i, j)] = Cell(i, j, target_aa, query_aa, diag,
                                                (i-1, j-1), score_diag)
                elif direct == up:
                    self.dpMat[(i, j)] = Cell(i, j, target_aa, "", up,
                                                (i-1, j), score_up)
                elif direct == left:
                    self.dpMat[(i, j)] = Cell(i, j, "", query_aa, left,
                                                (i, j-1), score_left)
    def fill_local(self):
        ncols = self.ncols #len(self.query)
        nrows = self.nrows #len(self.target)
        restart, diag, up, left = range(4)

        self.dpMat[(0, 0)] = Cell(0,0,None, None, None, None, 0.0)
        for j in range(1, ncols):
            # pad cols
            query_aa = self.query[j-1].upper()
            val = 0.0
            self.dpMat[(0, j)] = Cell(0,j,"", query_aa, left, (0,j-1), val)
            #self.dpMat[(0, j)] = Cell(0,j,None, None, None, None, val)
        for i in range(1,nrows):
            # pad rows
            val = 0.0
            target_aa = self.target[i-1].upper()
            self.dpMat[(i, 0)] = Cell(i,0, target_aa, "", up, (i-1, 0), val)
            #self.dpMat[(i, 0)] = Cell(i,0, None, None, None, None, val)

        #Fill the score matrix
        for i in range(1, nrows):
            target_aa = self.target[i-1].upper()
            for j in range(1, ncols):
                query_aa = self.query[j-1].upper()
                # no gap: diag
                score_diag = self.dpMat[(i-1, j-1)].value + \
                             self.sub_mat.get_value(target_aa, query_aa)
                # up
                score_up = self.dpMat[(i-1, j)].value + \
                            self.affine_gap_penalty(i-1)
                # left
                score_left = self.dpMat[(i, j-1)].value + \
                            self.affine_gap_penalty(j-1)

                max_score, direct = self._argmax([restart, score_diag, 
                                                score_up, score_left])

                if max_score > self.best_score:
                    self.best_score = max_score
                    self.best_score_inx = (i,j)

                if direct == diag:
                    self.dpMat[(i, j)] = Cell(i, j, target_aa, query_aa, diag,
                                                (i-1, j-1), score_diag)
                elif direct == up:
                    self.dpMat[(i, j)] = Cell(i, j, target_aa, "", up,
                                                (i-1, j), score_up)
                elif direct == left:
                    self.dpMat[(i, j)] = Cell(i, j, "", query_aa, left,
                                                (i, j-1), score_left)
                elif direct == restart:
                    self.dpMat[(i, j)] = Cell(i, j, "", "", left,
                                                None, 0)
    def get_cell(self, row, col):
        if row >= self.nrows:
            print("Error: row is overflow({})".format(self.nrows))
            sys.exit(1)
        if col >= self.ncols:
            print("Error: col is overflow({})".format(self.ncols))
            sys.exit(1)
        return self.dpMat[(row, col)]

    def affine_gap_penalty(self, length, penalize_extend=True):
        if length <=0:
            return 0
        # Affine Gap Function: gamma(n) = -1*(d + (n-1)*e)
        # where, d = gap_open and e = gap_extension
        # 
        pe = self.gap_open + self.gap_extension * length
        if not penalize_extend:
            pe -= self.gap_extension
        return -1*pe

    def print_dpMat(self):
        """
            Print dp matrix
        Args:

        Kwargs:
            header (list or str): 
        """
        for i in range(self.nrows):
            r = ""
            for j in range(self.ncols):
                val = self.dpMat[(i,j)].value
                r += "{:8.3f} ".format(val)
            print(r)

class Alignment(object):
    """
        Base class for pairwised sequence alignment
    .. note::
        This base class uses two patterns, composite and registry
    """

    __metaclass__ = ClassRegistry

    def __init__(self, *args, **kwargs):
        """Initialization.

        Args:
           algo (class):

        Kwargs:
        """
        super(Alignment, self).__init__()
        algo = kwargs.pop('algo', None)
        if algo is None:
            return
        #if algo:
            # instance, developped by a user or subclasses of Alignment
        #    self.__algo = algo()
        #elif self.__isstr(algo):
            # is string
        #    self.__algo = self.__create(algo, *args, **kwargs)
        self.__algo = self.__create(algo, *args, **kwargs)

    def __create(self, clsname, *args, **kwargs):
        """ Create an instance for a given class in str or name

        """
        obj = None
        for cls in self.registry:
            if clsname == cls.__name__ or clsname == cls:
                obj = cls(*args, **kwargs)
                break
        if not obj:
            print("Unknown class {}".format())
            sys.exit(1)
        #if obj:
            # Initialize object
        #    obj.__init__(*args, **kwargs)
        #else:
        #    print("Unknown class {}".format())
        return obj

    def __str__(self):
        #OK:    return self.__algo.__str__()
        #Failed return self.__getattr__(self.__algo, '__str__')
        algo_func_str_ = self.__getattr__(self.__algo, '__str__')
        return algo_func_str_()

    def __repr__(self):
        #return self.__algo.__repr__()
        algo_func_str_ = self.__getattr__(self.__algo, '__repr__')
        return algo_func_str_()

    def __call__(self, *args, **kwargs):
        return self.__algo(*args, **kwargs)

    def __getattr__(self, obj, attr):
        """ get delegation to the object """
        #return getattr(obj, attr)
        try:
            return self.__algo.__getattribute__(attr)
        except AttributeError:
            raise AttributeError('{0} object has no attribute `{1}`'
                .format(self.__class__.__name__, attr))

    def __isstr(self, s):
        try:
            return isinstance(s, basestring)
        except NameError:
            return isinstance(s, str)


class NWAlign(Alignment):
    """
        Implementation of Needleman-Wunsch algorithm
    """
    def __init__(self, *args, **kwargs):
        """Initialization.

        Args:
           param (type):

        Kwargs:

        """
        super(NWAlign, self).__init__(*args, **kwargs)
        
    def __str__(self):
        """ Serialize """
        return "Algo.name: {}".format(self.__class__.__name__)

    def __repr__(self):
        return self.__class__.__name__

    def __call__(self, *args, **kwargs):
        """
            Concrete implementation

        Args:
            query  (str): args[0]
            target (str): args[1]
            query_name (str): args[2], default is 'query' if not given
            target_name (str) : args[3], default is 'target' if not given

        Kwargs:
            sub_mat (biTK.SubMatrice.MatInfo.SubstitutionMatrix): default is NUC42
            gap_open (positive int): default is 10.0
            gap_extension (positive int): default is 1.0
        """
        nargc = len(args)
        if nargc == 4:
            query, target, query_name, target_name = args
        elif nargc == 3:
            query, target, query_name = args
            target_name = 'target'
        elif nargc == 2:
            query, target = args
            query_name = 'query'
            target_name = 'target'
        else:
            print("Error: wrong number of input")
            sys.exit(1)

        gap_open = kwargs.pop('gap_open', 10.0)
        gap_extension = kwargs.pop('gap_extension', 1.0)
        gap_symbol = kwargs.pop('gap_symbol','-')
        align_type = kwargs.pop('align_type', 'nt')
        if not align_type in ['nt', 'pr']:
            print("Error: unknown specified align type '{}', (nt or pr)."
                  .format(align_type))
            sys.exit(1)
        subMat_name = kwargs.pop('substitutionMatrix', 'nuc42')
        subMat = SUBMat(subMat_name)
        if not isinstance(subMat, SUBMat):
            print("Error: input type is not biTK.ngs.SubMatrice.MatInfo.SubstitutionMatrix")
            sys.exit(1)
        if not query or not target:
            return []
        
        # fill the DP matrix and traceback
        dpMatrix = DynamicProgMatrix(subMat, target, query, align_type, 
                        gap_open, gap_extension, gap_symbol=gap_symbol)
        # fill the matrix
        dpMatrix.fill(strategy='global')

        # trace back from the bottom right corner
        aligned_target = ""
        aligned_query = ""
        current_ = dpMatrix.get_cell(len(target), len(query))
        next_ = dpMatrix.get_cell(current_.parent[0], current_.parent[1])
        # by definition, the bottom corner cell is the optimal score 
        # for the alignment
        aligned_score = current_.value

        while next_:
            if ((next_.col == current_.col - 1 ) and 
                        (next_.row == current_.row - 1)):
                #1. move in the diagnal direction
                aligned_target += current_.target_aa
                aligned_query  += current_.query_aa
            elif ((next_.col == current_.col) and 
                        (next_.row == current_.row - 1)):
                #2. move in the upward direction
                aligned_target += current_.target_aa
                aligned_query  += gap_symbol
            elif ((next_.col == current_.col - 1) and 
                        (next_.row == current_.row)):
                #3. move in the left direction
                aligned_target += gap_symbol
                aligned_query  += current_.query_aa
            # move to next level
            current_ = next_
            if current_.parent is None:
                next_ = None
            else:
                next_ = dpMatrix.get_cell(current_.parent[0], current_.parent[1])
        now = datetime.datetime.now()
        # reverse and save to alignInfo
        alignInfo = AlignInfo(query_name, target_name,
                            query, target, 
                            aligned_score, 
                            aligned_query[::-1], 
                            aligned_target[::-1],
                    method='Global', subMatrix=subMat.name,tbMatrix=dpMatrix,
                    gap_open=gap_open, gap_extension=gap_extension,
                    gap_symbol=gap_symbol, align_type=align_type, 
                    t="{0:%Y-%m-%d %H:%M:%S}".format(now))

        return alignInfo

class SWAlign(Alignment):
    """
        Implementation of Smith-Watermann algorithm
    """
    def __init__(self, *args, **kwargs):
        """Initialization.

        Args:
           param (type):

        Kwargs:

        """
        super(SWAlign, self).__init__(*args, **kwargs)

    def __str__(self):
        """ Serialize """
        return "Algo.name: {}".format(self.__class__.__name__)

    def __repr__(self):
        return self.__class__.__name__

    def __call__(self, *args, **kwargs):
        """
            Concrete implementation

        Args:
            query  (str): args[0]
            target (str): args[1]
            query_name (str): args[2], default is 'query' if not given
            target_name (str) : args[3], default is 'target' if not given

        Kwargs:
            sub_mat (biTK.SubMatrice.MatInfo.SubstitutionMatrix): default is NUC42
            gap_open (positive int): default is 10.0
            gap_extension (positive int): default is 1.0
        """
        nargc = len(args)
        if nargc == 4:
            query, target, query_name, target_name = args
        elif nargc == 3:
            query, target, query_name = args
            target_name = 'target'
        elif nargc == 2:
            query, target = args
            query_name = 'query'
            target_name = 'target'
        else:
            print("Error: wrong number of input")
            sys.exit(1)

        gap_open = kwargs.pop('gap_open', 10.0)
        gap_extension = kwargs.pop('gap_extension', 1.0)
        gap_symbol = kwargs.pop('gap_symbol','-')
        align_type = kwargs.pop('align_type', 'nt')
        if not align_type in ['nt', 'pr']:
            print("Error: unknown specified align type '{}', (nt or pr)."
                  .format(align_type))
            sys.exit(1)
        subMat_name = kwargs.pop('substitutionMatrix', 'nuc42')
        subMat = SUBMat(subMat_name)
        if not isinstance(subMat, SUBMat):
            print("Error: input type is not biTK.ngs.SubMatrice.MatInfo.SubstitutionMatrix")
            sys.exit(1)
        if not query or not target:
            return []
        
        # fill the DP matrix and traceback
        dpMatrix = DynamicProgMatrix(subMat, target, query, align_type, 
                        gap_open, gap_extension, gap_symbol=gap_symbol)
        # fill the matrix
        dpMatrix.fill(strategy='local')

        # trace back from the bottom right corner
        aligned_target = ""
        aligned_query = ""
        best_score_inx = dpMatrix.best_score_inx
        current_ = dpMatrix.get_cell(best_score_inx[0], best_score_inx[1])
        next_ = dpMatrix.get_cell(current_.parent[0], current_.parent[1])
        # by definition, the maximum score is the optimal score 
        # for the alignment
        aligned_score = dpMatrix.best_score

        while next_:
            if ((next_.col == current_.col - 1 ) and 
                        (next_.row == current_.row - 1)):
                #1. move in the diagnal direction
                aligned_target += current_.target_aa
                aligned_query  += current_.query_aa
            elif ((next_.col == current_.col) and 
                        (next_.row == current_.row - 1)):
                #2. move in the upward direction
                aligned_target += current_.target_aa
                aligned_query  += gap_symbol
            elif ((next_.col == current_.col - 1) and 
                        (next_.row == current_.row)):
                #3. move in the left direction
                aligned_target += gap_symbol
                aligned_query  += current_.query_aa
            # move to next level
            current_ = next_
            if current_.parent is None:
                next_ = None
            else:
                next_ = dpMatrix.get_cell(current_.parent[0], current_.parent[1])
        now = datetime.datetime.now()

        # reverse and save to alignInfo
        alignInfo = AlignInfo(query_name, target_name,
                            query, target, 
                            aligned_score, 
                            aligned_query[::-1], 
                            aligned_target[::-1],
                    method='Local', subMatrix=subMat.name, tbMatrix=dpMatrix,
                    gap_open=gap_open, gap_extension=gap_extension,
                    gap_symbol=gap_symbol, align_type=align_type, 
                    t="{0:%Y-%m-%d %H:%M:%S}".format(now))

        return alignInfo
