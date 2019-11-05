import sys

sys.path.append(".")

import Bpowspec
from pylab import *

close('all')

size3d = array([1.,1.,1.])

N = array([128,128,128],dtype=int32)

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
Pk = exp(-k**2/2/(kmax/5)**2)*cos(k/kmax*15)**2


plot(k,Pk)
xlabel('k')
ylabel('P(k)')

seed = 1233453


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


Bharmxobs = Bpowspec.map2harm(Bx,Delta)
Bharmyobs = Bpowspec.map2harm(By,Delta)
Bharmzobs = Bpowspec.map2harm(Bz,Delta)


kobsbin = linspace(0,kmax,50)
Pkobs = Bpowspec.harm2Pk(Bharmxobs,Bharmyobs,Bharmzobs, Deltak, kobsbin)
Pkobs2 = Bpowspec.harm2Pk(Bharmx,Bharmy,Bharmz, Deltak, kobsbin)

kobs = kobsbin[:-1] + diff(kobsbin)/2.

figure(1)
step(kobs,Pkobs,where='mid',label='though map')
step(kobs,Pkobs2,where='mid',label='direct harm')
legend()

figure(2)
imshow(Bz[0])

figure(3)
imshow( real(Bharmx - Bharmxobs)[0])                                    



show()

