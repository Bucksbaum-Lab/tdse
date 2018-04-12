# Python [matplotlib] animation technique with courtesy to Dr. Vanderplas
#
import matplotlib
matplotlib.use('TKAgg')
from numpy import pi, sqrt, exp, tanh
import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
from schr2d import schr2d
import Conversion as conv

def gauss(x, y, ax, ay, x0, y0, phase):
    """
    2D gaussian wave packet
 
    psi_0(x) = N_0 x exp[-0.5 a^2 x^2] x H_0(x)
             = (a / sqrt(pi))**0.5 x exp [-0.5 a^2 x^2]
    a = sqrt(mx wx / hbar)
    """
    g_x = ((ax / sqrt(pi)) ** 0.5
          * exp(-0.5 * ((x - x0) * ax) ** 2))
    g_y = ((ay / sqrt(pi)) ** 0.5
          * exp(-0.5 * ((y - y0) * ay) ** 2))

    gxy = np.zeros((len(x),len(y)), dtype=float)
    for i, _gx in enumerate(g_x):
        for j, _gy in enumerate(g_y):
            gxy[i,j] = _gx * _gy           

    gxy2 = (1.0 / sqrt(1.0+abs(phase))) * np.array([gxy, phase*gxy], dtype=float) 

    return gxy2

######################################################################
#
# 2D CoIn model [au]
#
######################################################################
dt = 5.0
t_max = 12000.0
frames = int(t_max / dt)

hbar = 1.0 
mx   = 1.0 * conv.am2au    
my   = 5.0 * conv.am2au
kap  = 0.05
lmb  = 0.2
x0   = 2.5 # [au]
y0   = 0.0
wx   = 3000.0 * conv.ic2au
wy   = 2000.0 * conv.ic2au
#
alp_x = sqrt(mx * wx / hbar) 
alp_y = sqrt(my * wy / hbar)

# prepare calculation
N  = 2**7
dx = 16.0 / float(N)
dy = dx
x  = dx * (np.arange(N) - 0.5 * N)
y  = dy * (np.arange(N) - 0.5 * N)

# diabatic potential
Vp  = np.zeros((N,N), dtype=float)
Vm  = np.zeros((N,N), dtype=float)
V12 = np.zeros((N,N), dtype=float)
for i, _x in enumerate(x):
    for j, _y in enumerate(y):
        Vp[i][j]  = -kap * (_x - x0) * 0.5 * (1.0 - tanh(_x-x0)) + 0.5 * my * wy**2 * _y**2
        Vm[i][j]  =  kap * (_x + x0) * 0.5 * (1.0 + tanh(_x+x0)) + 0.5 * my * wy**2 * _y**2
        V12[i][j] =  lmb * _y

V_xy = [Vp, Vm, V12]

# for potential plot
Vp_x = np.zeros(N, dtype=float)
Vm_x = np.zeros(N, dtype=float)
V_y  = np.zeros(N, dtype=float)
for i in range(N):
    Vp_x[i] = Vp[i][int(N/2)] * conv.au2ev
    Vm_x[i] = Vm[i][int(N/2)] * conv.au2ev
    V_y[i]  = 0.5 * my * wy**2 * y[i]**2 * conv.au2ev

# initial wave packet
phase = 0 #-1.0
  
psi0 = gauss(x, y, alp_x, alp_y, 0.0, 0.0, phase)
 
S = schr2d (x=x,
            y=y,
            psi0=psi0,
            V_xy=V_xy,
            mx=mx,
            my=my,
            nstates=2)

######################################################################
#
# Construct plot
#
######################################################################
fig = plt.figure()

# 
xlim = (-8, 8)
ylim = (-3, 3)

# top axes show the x-space data
ymin = 0.0
ymax = 7.0 
ax1 = fig.add_subplot(211, xlim=xlim,
                      ylim=(ymin - 0.2 * (ymax - ymin),
                            ymax + 0.2 * (ymax - ymin)))
                          
psi_x_line, = ax1.plot([], [], c='r', label=r'$|\psi(x)|$')
Vp_x_line,  = ax1.plot([], [], c='k', label=r'$V_+(x)$')
Vm_x_line,  = ax1.plot([], [], c='b', label=r'$V_-(x)$')

title = ax1.set_title("")
ax1.legend(prop=dict(size=12))
ax1.set_xlabel('$x$')
ax1.set_ylabel(r'$|\psi(x)|$')

# bottom axes show the y-space data
ymin = 0.0 
ymax = 3. 
ax2 = fig.add_subplot(212, xlim=ylim,
                      ylim=(ymin - 0.2 * (ymax - ymin),
                            ymax + 0.2 * (ymax - ymin)))
psi_y_line, = ax2.plot([], [], c='r', label=r'$|\psi(y)|$')
V_y_line,   = ax2.plot([], [], c='k', label=r'$V(y)$')

ax2.legend(prop=dict(size=12))
ax2.set_xlabel('$y$')
ax2.set_ylabel(r'$|\psi(y)|$')

######################################################################
#
# Animation
#
######################################################################

def init():
    psi_x_line.set_data([], [])
    psi_y_line.set_data([], [])
    title.set_text("")
    return (psi_x_line, psi_y_line, title)


def animate(i):
    S._step(dt)
    # reduced rho(x)
    psi_x_line.set_data(S.x, 6.0*abs(S.psi_x))
    # V(x)
    Vp_x_line.set_data(S.x, Vp_x)
    Vm_x_line.set_data(S.x, Vm_x)
    # reduced rho(y)
    psi_y_line.set_data(S.y, abs(S.psi_y))
    # V(y)
    V_y_line.set_data(S.y, V_y)
    title.set_text("t = %.2f au" % S.t)
    return (psi_x_line, psi_y_line, title)

# Animation
# For Linux 
#
# blit=True - only re-draw the changed parts.
#
#anim = animation.FuncAnimation(fig, animate, init_func=init,
#                               frames=frames, interval=30, blit=True)
# For MacOS
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=frames, interval=30, blit=False)



# uncomment the following line to save the video in mp4 format.  This
# requires either mencoder or ffmpeg to be installed on your system
#anim.save('2D_CoIn.mp4', fps=15,
#          extra_args=['-vcodec', 'libx264'])

plt.show()
