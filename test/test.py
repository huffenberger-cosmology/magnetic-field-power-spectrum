import sys

sys.path.append(".")

import Bpowspec
from pylab import *

close('all')

size3d = array([1.,1.,1.])

N = array([256,256,256],dtype=int32)

Delta = size3d/N

Deltak = Bpowspec.Delta2k(Delta,N)

print(Delta)
print(Deltak)


kvec = array([1.,0.2,0.5])
khat = kvec/norm(kvec)

P = Bpowspec.build_projector(khat)


print(khat)
print(P)


kmax = amax(N*Deltak)
k = linspace(0,kmax, 500);
Pk = exp(-k**2/2/(kmax/40)**2)


plot(k,Pk)
xlabel('k')
ylabel('P(k)')

seed = 1233


kx,ky,kz = Bpowspec.kvecs(N,Deltak)

Bharmx,Bharmy,Bharmz = Bpowspec.Pk2harm(k,Pk,N,kmax,Deltak, seed)

DivBharm = kx*Bharmx + ky*Bharmy + kz*Bharmz


print("==============\nDeltak * max(norm(DivBharm)) = %e" % (Deltak[0] * amax(norm(DivBharm))))
print("rms(norm(Bharm)) = %f\n==============" % std(sqrt(real(
    Bharmx*conj(Bharmx) +
    Bharmy*conj(Bharmy) +
    Bharmz*conj(Bharmz) ))))



Bx = Bpowspec.harm2map(Bharmx,Delta)
By = Bpowspec.harm2map(Bharmy,Delta)
Bz = Bpowspec.harm2map(Bharmz,Delta)


figure()
imshow(Bz[0])
show()

