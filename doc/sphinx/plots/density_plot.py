import matplotlib.pyplot as plt
import periodictable

params = {'figure.figsize': (10,8)}
plt.rcParams.update(params)

D = [(el.number,el.density,el.symbol)
    for el in periodictable.elements]

bbox = dict(boxstyle="round",lw=1,ec=(0,0,0),fc=(0.85,0.8,0.8))
for Z,density,sym in D:
    if density is not None:
        plt.text(Z,density,sym,bbox=bbox)
plt.axis([0,110,0,25])
plt.xlabel('Element number')
plt.ylabel('Density of element')
plt.title('Density for elements')
