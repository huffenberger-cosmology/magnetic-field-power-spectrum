import sys

sys.path.append(".")

import Bpowspec
from pylab import *

size3d = array([1.,1.,1.])

N = array([128,128,128],dtype=int32)

Delta = size3d/N

Deltak = Bpowspec.Delta2k(Delta,N)

print(Delta)
print(Deltak)
