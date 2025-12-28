
import numpy as np
from suspensionbridge import *

N = 20

sb = SuspensionBridge_xz(N)
#sb.prune_nu() # remove null mode
print(sb.Ev)

sb.plot_spectrum("spectrum_XZ_N{}.pdf".format(N))
for mode in range(8):
 sb.plot_XZ_modes(mode, "lowest_N{}_XZ_modes{}.pdf".format(N,mode))
print(sb.nu)       

L=0
for i in range(1,N):
  L+=np.cos(sb.thetas[i])
print(L)

h=0
for i in range(N//2,N):
  h+=np.sin(sb.thetas[i])
print(h)
