# Wave packet interferometry on top of a CoIn
import sys
import numpy as np
from numpy import pi, exp, sqrt, cos, sin, conjugate
from scipy.fftpack import fft, ifft
from numpy import linalg

class schr2d(object):
    def __init__(self, x, y, psi0, V_xy, mx, my, nstates):
        self.x, self.y, self.psi0, self.V_xy = map(np.array, (x, y, psi0, V_xy))
        self.Nx = len(self.x)
        self.Ny = len(self.y) 

        self.mx  = mx
        self.my  = my
        self.nstates = nstates
        self.t   = 0.0
        
        self.dx  = x[1] - x[0]
        self.dy  = y[1] - y[0]   
        self.dkx = 2. * pi / (self.Nx * self.dx)
        self.dky = 2. * pi / (self.Ny * self.dy)

        self.k0x = -0.5 * self.Nx * self.dkx
        self.k0y = -0.5 * self.Ny * self.dky

        self.kx = self.k0x + self.dkx * np.arange(self.Nx)
        self.ky = self.k0y + self.dky * np.arange(self.Ny)

        # psi0[nstates][Nx][Ny]
        self.psi0  = np.array(psi0, dtype=complex)
        self.psi   = self.psi0.copy()

        self.wf_norm_x_y(self.psi)

        #
        self.psi_x = None
        self.psi_y = None
         
    def wf_norm_x_y(self, wavefunc):
        """
        compute reduced rho(x) and rho(y) 
                        self.psi_x and self.psi_y
        """
        Nx = len(self.x)
        Ny = len(self.y)
        psi_x = np.zeros(Nx,dtype=float)
        psi_y = np.zeros(Ny,dtype=float)

        psi_xy = np.zeros((Nx,Ny), dtype=float)
        norm = 0.0
        for i in range(Nx):
            for j in range(Ny):
                for s in range(self.nstates):
                    psi_xy[i,j] += wavefunc[s,i,j] * conjugate(wavefunc[s,i,j])
                    norm        += wavefunc[s,i,j] * conjugate(wavefunc[s,i,j])
                #
                psi_xy[i,j] = sqrt(psi_xy[i,j])
                #
        norm = norm * self.dx * self.dy
        #sys.stdout.write("norm = %.3f\n" % (norm))
        print "norm=", norm
        #
        # compute rho(x)
        for i in range(Nx):
            for j in range(Ny):
                psi_x[i] += psi_xy[i,j]
            psi_x[i] = psi_x[i] * self.dx
        #
        # compute rho(y)
        for j in range(Ny):
            for i in range(Nx):
                psi_y[j] += psi_xy[i,j]
            psi_y[j] = psi_y[j] * self.dy
    
        self.psi_x = psi_x.copy()
        self.psi_y = psi_y.copy()


    def _step(self, dt):
        #
        # first with potential operator
        self.potentOp(dt)

        # second with kinetical operator
        self.kinOp(dt)

        # another part of potential operator 
        self.potentOp(dt)

        self.t += dt
        self.wf_norm_x_y(self.psi)


    def kinOp(self, dt):
        psi = self.psi.copy()
        Nx  = self.Nx
        Ny  = self.Ny
        for s in range(self.nstates):
            #
            # exp[-iTdt][psi] for x
            #
            for iy in range(Ny):
                # forward FFT
                arr = np.zeros(Nx, dtype=complex)
                for ix in range(Nx): 
                    arr[ix] = psi[s, ix, iy] 

                arr = fft(arr)

                # exp[iTdt]
                for i, k in enumerate(self.kx):
                    arr[i] = arr[i] * exp( -1j * k**2 * dt / (2.0 * self.mx) ) 

                # backward FFT
                arr = ifft(arr)
               
                for ix in range(Nx):
                    psi[s,ix,iy] = arr[ix]

            #
            # exp[-iTdt][psi] of y
            #
            for ix in range(Nx):
                # forward FFT
                arr = np.zeros(Ny, dtype=complex) 
                for iy in range(Ny):
                    arr[iy] = psi[s,ix,iy]  
                
                arr = fft(arr)
 
                # exp[iTdt]
                for i, k in enumerate(self.ky):
                    arr[i] = arr[i] * exp( -1j * k**2 * dt / (2.0 * self.my) )

                # backward FFT
                arr = ifft(arr)
# 
                for iy in range(Ny):
                    psi[s,ix,iy] = arr[iy]
#
        self.psi = psi.copy()

    def potentOp(self, dt):
        """
        potential operator for dt/2
        """
        dth = 0.5 * dt 
        upot = np.zeros((self.nstates,self.nstates), dtype=complex)
        psi = self.psi.copy()
        #
        for ix in range(self.Nx):
            for iy in range(self.Ny):
                upot[0,0] = self.V_xy[0][ix][iy]
                upot[1,1] = self.V_xy[1][ix][iy]
                upot[0,1] = self.V_xy[2][ix][iy]
                upot[1,0] = upot[0,1]
                
                eig  = linalg.eig(upot)
                en   = exp(-1j * eig[0] * dth) # eigenwerte
                en   = np.diag(en)
                #print en
                upot = eig[1]                  # eigenvektoren

                cw1   = np.dot(en, upot)
                opmat = np.dot(upot.transpose(), cw1)  

                # ___chk___ unitarity
                #u = np.dot(opmat.conjugate(),opmat)
                #print u
                #exit(0)
                        
                # exp[-iVdt/2][psi]
                self.psi[0,ix,iy] = opmat[0,0] * psi[0,ix,iy] + opmat[0,1] * psi[1,ix,iy] 
                self.psi[1,ix,iy] = opmat[1,0] * psi[0,ix,iy] + opmat[1,1] * psi[1,ix,iy]

        #self.wf_norm_x_y(self.psi)


