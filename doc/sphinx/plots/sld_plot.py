import matplotlib.pyplot as plt
import periodictable.nsf

params = {'figure.figsize': (10,8)}
plt.rcParams.update(params)

periodictable.nsf.sld_plot()
