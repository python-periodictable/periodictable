import numpy as np
import matplotlib.pyplot as plt
import periodictable

fig_width = 10
fig_height = 8
fig_size =  [fig_width,fig_height]
params = {'figure.figsize': fig_size}
plt.rcParams.update(params)
Fe_2 = periodictable.Fe.ion[2]
Q = np.linspace(0,16,200)
M = Fe_2.magnetic_ff[Fe_2.charge].j0_Q(Q)
plt.xlabel(r'Magnetic Form Factor for Fe')
plt.ylabel(r'$\AA^{-1}$')
plt.title('Ion specific property for Fe')
plt.plot(Q,M)
