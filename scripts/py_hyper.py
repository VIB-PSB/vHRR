from ctypes import cdll, c_double
from os.path import dirname, abspath
lib = cdll.LoadLibrary(dirname(abspath(__file__))+"/hyper.so")

class Hyper(object):
    def __init__(self):
        lib.getHyperP.restype = c_double
        lib.getUpperCumulativeHyperP.restype = c_double
        lib.getLowerCumulativeHyperP.restype = c_double
        self.obj = lib.Hyper_new()
    
    def getHyperP(self, totalSize, sampleSize, totalHits, sampleHits):
        return lib.getHyperP(self.obj, totalSize, sampleSize, totalHits, sampleHits)
    
    def getUpperCumulativeHyperP(self, totalSize, sampleSize, totalHits, sampleHits):
        return lib.getUpperCumulativeHyperP(self.obj, totalSize, sampleSize, totalHits, sampleHits)
    
    def getLowerCumulativeHyperP(self, totalSize, sampleSize, totalHits, sampleHits):
        return lib.getLowerCumulativeHyperP(self.obj, totalSize, sampleSize, totalHits, sampleHits)

