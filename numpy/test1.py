#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

N = 32
Fs = 4.5e9
Fc = 0.12345e9

t = np.linspace(0,N-1,N)/Fs
s = np.sin(2*np.pi*Fc*t)

plt.plot(t, s, 'ro-')
plt.ylabel('value')
plt.xlabel('time')
plt.show()

#print('a = {}\n'.format(a))
#print(a.ndim, a.shape, a.size, a.dtype, a.itemsize, a)




