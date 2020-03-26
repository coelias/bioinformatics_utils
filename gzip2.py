import zlib 
from StringIO import StringIO
import sys

class gzfile:
    '''Module that implements a ultrafast gzip'''
    MAXSIZE=1024*64
    def __init__(self,path1=None,fileobj=None,sz=0):
        self.sz=sz
        assert not path1 or not fileobj
        self.dec = zlib.decompressobj(16+zlib.MAX_WBITS)
        if fileobj: self.fil1=fileobj
        else: self.fil1=open(path1,'rb')
        self.buff1=StringIO(self.dec.decompress(self.fil1.read(gzfile.MAXSIZE)))
        self._read=gzfile.MAXSIZE


    def __iter__(self):
        line=""
        while True:
            for line in self.buff1:
                if line[-1]!="\n":
                    break
                yield line

            if line and line[-1]=="\n": line=""
            infodec=self.fil1.read(gzfile.MAXSIZE)
            self._read+=gzfile.MAXSIZE
            if self.sz:
                sys.stderr.write('\r{0:.2f}% read'.format(100*(float(self._read)/self.sz)))

            info=line+self.dec.decompress(infodec)
            self.buff1=StringIO(info)
            if not infodec: break

        if line:
            yield line
    
    def seek(self,pos):
        self.dec = zlib.decompressobj(16+zlib.MAX_WBITS)
        self.fil1.seek(0)
        while True:
            infodec=self.fil1.read(gzfile.MAXSIZE)
            if not infodec: 
                self.buff1=StringIO()
                break
            infodec=self.dec.decompress(infodec)
            if len(infodec)<pos:
                pos-=len(infodec)
            else:
                self.buff1=StringIO(infodec[pos:])
                break


    def close(self):
        self.fil1.close()

if __name__=="__main__":
    import sys
    a=gzfile(sys.argv[1])
    for i in a:
        sys.stdout.write(i)

