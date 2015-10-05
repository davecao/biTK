# -*- coding: utf-8 -*-

__all__ = ['Buffer']

from gzip import GzipFile
from bz2 import BZ2File
from collections import deque


CHUNK = 16 * 1024


class Buffer(object):
    """ File buffer """
    def __init__ (self):
        super(Buffer, self).__init__()
        self.__buf = deque()
        self.__size = 0

    def __len__ (self):
        return self.__size

    def write(self, data):
        self.__buf.append(data)
        self.__size += len(data)

    def read(self, size=-1):
        if size < 0: size = self.__size
        ret_list = []
        while size > 0 and len(self.__buf):
            s = self.__buf.popleft()
            size -= len(s)
            ret_list.append(s)
        if size < 0:
            ret_list[-1], remainder = ret_list[-1][:size], ret_list[-1][size:]
            self.__buf.appendleft(remainder)
        ret = ''.join(ret_list)
        self.__size -= len(ret)
        return ret

    def flush (self):
        pass

    def close (self):
        pass

class Bz2CompressReadStream(object):
    def __init__ (self, fileobj):
        self.__input = fileobj
        self.__buf = Buffer()
        self.__bz2 = BZ2File(None, mode='rb', fileobj=self.__buf)

    def read(self, size=-1):
        while size < 0 or len(self.__buf) < size:
            s = self.__input.read(CHUNK)
            if not s:
                self.__gzip.close()
                break
            self.__gzip.write(s)
        return self.__buf.read(size)

class GzipCompressReadStream (object):

    def __init__ (self, fileobj):
        self.__input = fileobj
        self.__buf = Buffer()
        self.__gzip = GzipFile(None, mode='wb', fileobj=self.__buf)

    def read(self, size=-1):
        while size < 0 or len(self.__buf) < size:
            s = self.__input.read(CHUNK)
            if not s:
                self.__gzip.close()
                break
            self.__gzip.write(s)
        return self.__buf.read(size)

class MemoryGzipFile(gzip.GzipFile):
    """
    A GzipFile subclass designed to be used with in memory file like
    objects, i.e. StringIO objects.
    """
    def write_crc_and_filesize(self):
        """
        Flush and write the CRC and filesize. Normally this is done
        in the close() method. However, for in memory file objects,
        doing this in close() is too late.
        """
        self.fileobj.write(self.compress.flush())
        gzip.write32u(self.fileobj, self.crc)
        # self.size may exceed 2GB, or even 4GB
        gzip.write32u(self.fileobj, self.size & 0xffffffffL)

    def close(self):
        if self.fileobj is None:
            return
        self.fileobj = None
        if self.myfileobj:
            self.myfileobj.close()
            self.myfileobj = None
