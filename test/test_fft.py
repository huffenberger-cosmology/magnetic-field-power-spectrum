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


kmax = amax(N*Deltak)
k = linspace(0,kmax, 500);
Pk = exp(-k**2/2/(kmax/5)**2)*cos(k/kmax*15)**2


figure(1)
plot(k,Pk)
xlabel('k')
ylabel('P(k)')

seed = 1233

harm = Bpowspec.scalarPk2harm(k,Pk,N,kmax,Deltak, seed)

m = Bpowspec.harm2map(harm,Delta)

harmobs = Bpowspec.map2harm(m,Delta)

figure(2)
imshow(m[0])

figure(3)
imshow( real(harm - harmobs)[0])                                    
colorbar()

print("rms real diff = %e" % std(real(harm - harmobs)))
print("rms imag diff = %e" % std(imag(harm - harmobs)))

figure(4)
plot(harm[0,:,0], label='harm')
plot(harmobs[0,:,0], label='harmobs')
legend()
     

show()

