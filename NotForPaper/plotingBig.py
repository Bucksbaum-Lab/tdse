import scipy.io
import matplotlib as mpl
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

print 'importing file'

tdsedat = scipy.io.loadmat('C:\Chelsea\TDSE\OutputOct1-2015.mat')

print 'making figure'

fig = plt.figure()
ax=fig.gca(projection='3d')
x1=tdsedat['x1']
x2=tdsedat['x2']
wavefunctionE = np.abs(tdsedat['wavefunctionE'])
wavefunctionG = np.abs(tdsedat['wavefunctionG'])

print 'ploting figure'

ax.plot_surface(x1,x2,wavefunctionE, cmap=cm.coolwarm, linewidth = 0)
plt.show()
