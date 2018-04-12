##  Copyright (C) 2010 O. Vendrell
##  QDTK, a set of modules for quantum dynamics and methods development
##
##  This file is part of QDTK.
##
##  QDTK is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
## $Id: Conversion.py 1220 2014-09-02 21:25:23Z ovendrel $

"""
Compendium of conversion factors and physical constants:

PI                = 3.1415926535897    #Pi
E                 = 2.718281828459     #E
FPC_c             = 299792458          #speed of light in vacuum (def) m/s
FPC_h_BAR         = 1.05457266e-34     #Planck constant, reduced (63) J s
FPC_h_BAR_MeVs    = 6.582122e-22       #Planck constant, reduced (20) MeV s
FPC_e_C           = 1.60217733e-19     #electron charge magnitude (49) C
FPC_e_ESU         = 4.8032068e-10      #electron charge magnitude (15) esu
FPC_hBARc         = 197.327053         #conversion constant hbar*c (59) MeV Fm
FPC_hBARc2        = 0.38937966         #conversion constant (hbar*c)^2 (23) GeV^2 mbarn
FPC_m_e_kg        = 9.1093897e-31      #electron mass (54) kg
FPC_m_e_MeV       = 0.51099906         #electron mass (15) MeV/c^2
FPC_m_P_MeV       = 938.27231          #proton mass (28) MeV/c^2
FPC_m_P_u         = 1.00727647         #proton mass (12) u
FPC_m_P_kg        = 1.6726231e-27      #proton mass (10) kg
FPC_m_P_M_E       = 1836.152701        #proton mass (37) m_e
FPC_m_D_MeV       = 1875.61339         #deuteron mass (57) MeV/c^2
FPC_u_MeV         = 931.49432          #unified atomic mass unit (u) (28) MeV/c^2
FPC_u_kg          = 1.6605402e-27      #unified atomic mass unit (u) (10) kg
FPC_EPSILON_0     = 8.854187817e-12    #permittivity of free space F/m
FPC_MU_0          = 1.2566370614e-06   #permeability of free space N/A^2
FPC_ALPHA         = 0.0072973530796448 #fine-structure constant (61)
FPC_r_e           = 2.81794092e-15     #classical electron radius (38) m
FPC_LAMBDA_BAR_e  = 3.86159323e-13     #electron Compton wavelength (35) m
FPC_a_0           = 5.29177249e-11     #Bohr radius (mnucleus= infty) (24) m
FPC_LAMBDA_1EV    = 1.23984244e-06     #wavelength of 1 eV/c particle (37) m
FPC_R_INFINITY_EV = 13.6056981         #Rydberg energy (mnucleus = infinity) (40) eV
FPC_SIGMA_0_BARN  = 0.66524616         #Thomson cross section (18) barn
FPC_MU_B_MeV_T    = 5.78838263e-11     #Bohr magneton (52)  MeV/T
FPC_MU_N_MeV_T    = 3.15245166e-14     #nuclear magneton (28) MeV/T
FPC_E_M_e         = 175881962000       #electron cyclotron freq./field (53) C/kg (rad/sT)
FPC_E_M_P         = 95788309           #proton cyclotron freq./field (29) C/kg (rad/sT)
FPC_G_SI          = 6.67259e-11        #gravitational constant (85) m^3/kgs^2
FPC_G_P           = 6.70711e-39        #gravitational constant (86) h_bar c (GeV/c^2)^{-2}
FPC_g             = 9.80665            #standard grav. accel., sea level m/s^2
FPC_N_A           = 6.0221367e+23      #Avogadro constant (36)  /mole
FPC_K_B           = 1.380658e-23       #Boltzmann constant (12) J/K
FPC_K_B_EV        = 8.617385e-05       #Boltzmann constant (73) eV/K
FPC_V_MOLAR       = 0.0224141          #molar volume, ideal gas at STP (19) m^3/mole
FPC_LAMBDAT       = 0.002897756        #Wien displacement law constant (24) m K
FPC_SIGMA_SB      = 5.67051e-08        #Stefan-Boltzmann constant (19) W/m^2K^4
FPC_G_F           = 1.16639e-05        #Fermi coupling constant (2)  GeV^{-2}
FPC_SIN2_THETA_W  = 0.23192            #weak mixing angle  (5)  at M_Z
FPC_M_W           = 80.22              #W boson mass (26) GeV/c^2
FPC_M_Z0          = 91.187             #Z_0 boson mass (7) GeV/c^2
FPC_G_S           = 0.117              #strong coupling constant (5))  at M_Z
AUCHARGE          = 1.60217733e-19     #1 a.u. charge in SI, elementary charge
AUMASS            = 9.1093897e-31      #1 a.u. mass in SI, electron rest mass
AUENERGY          = 4.3597482e-18      #1 a.u. energy in SI, hartree energy
AULENGTH          = 5.29177249e-11     #1 a.u. length in SI, bohrradius
AUACTION          = 1.05457266e-34     #1 a.u. action in SI, Planck constant/(2*PI)
AUTIME            = 2.4188843341e-17   #1 a.u. time in SI
AUFORCE           = 8.2387295e-08      #1 a.u. force in  SI
AUVELOCITY        = 2187691.42         #1 a.u. velocity in SI
AUMOMENTUM        = 1.9928534e-34      #1 a.u. momentum in SI
AUEFIELD          = 514220820000       #1 a.u. el. field in SI
AUEDIPOLE         = 8.4783579e-30      #1 a.u. el. dipole in SI
AUMFLUX           = 235051.808         #1 a.u. magn. flux in SI
AUMDIPOLE         = 1.85480308e-23     #1 a.u. magn. flux in  SI
ELECTRONVOLT      = 1.60217733e-19     #1 eV charge in SI
ATOMICMASS        = 1.6605402e-27      #1 amu mass in SI
ASTRONOMICALUNIT  = 149597870000       #1 AE length in SI
MOL               = 6.0221367e+23      #Avogadro constant
au2an             = 0.529177249        #[a.u.] -> [Angstroem]
an2au             = 1.8897259885789    #[Angstroem] -> [a.u.]
au2ev             = 27.211396131788    #[a.u.] -> [eV]
ev2au             = 0.036749308824762  #[eV] -> [a.u.]
au2am             = 0.00054857989586762#[a.u.] -> [atomar mass units]
am2au             = 1822.8885300626    #[atomar mass units] -> [a.u.]
au2fs             = 0.024188843341     #[a.u.] -> [fs]
fs2au             = 41.341373206755    #[fs] -> [a.u.]
ev2kc             = 23.05              #[eV] -> [kcal/mol]
kc2ev             = 0.043364476254697  #[kcal/mol] -> [eV]
au2kc             = 627.50431878767    #[a.u.] -> [kcal/mol]
kc2au             = 0.0015936145299079 #[kcal/mol] -> [a.u.]
au2ic             = 219474.7           #[au] -> [1/cm]
ic2au             = 4.5563338279993e-06#[1/cm] -> [au]
ev2ic             = 8065.54353         #[eV] -> [1/cm]
ic2ev             = 0.0001239842047    #[1/cm] -> [eV]
ev2kj             = 96.485             #[eV] -> [kJoule/mol]
kj2ev             = 0.0103643          #[kJoule/mol] -> [eV]
in2cm             = 2.54               #[inch] -> [cm]
ft2m              = 0.3048             #[feet] -> [cm]
yd2m              = 0.9144             #[yard] -> [cm]
cm2in             = 0.39370078740157   #[inch] -> [cm]
m2ft              = 3.2808398950131    #[feet] -> [cm]
m2yd              = 1.0936132983377    #[yard] -> [cm]
"""

PI                = 3.1415926535897    #Pi
E                 = 2.718281828459     #E
FPC_c             = 299792458          #speed of light in vacuum (def) m/s
FPC_h_BAR         = 1.05457266e-34     #Planck constant, reduced (63) J s
FPC_h_BAR_MeVs    = 6.582122e-22       #Planck constant, reduced (20) MeV s
FPC_e_C           = 1.60217733e-19     #electron charge magnitude (49) C
FPC_e_ESU         = 4.8032068e-10      #electron charge magnitude (15) esu
FPC_hBARc         = 197.327053         #conversion constant hbar*c (59) MeV Fm
FPC_hBARc2        = 0.38937966         #conversion constant (hbar*c)^2 (23) GeV^2 mbarn
FPC_m_e_kg        = 9.1093897e-31      #electron mass (54) kg
FPC_m_e_MeV       = 0.51099906         #electron mass (15) MeV/c^2
FPC_m_P_MeV       = 938.27231          #proton mass (28) MeV/c^2
FPC_m_P_u         = 1.00727647         #proton mass (12) u
FPC_m_P_kg        = 1.6726231e-27      #proton mass (10) kg
FPC_m_P_M_E       = 1836.152701        #proton mass (37) m_e
FPC_m_D_MeV       = 1875.61339         #deuteron mass (57) MeV/c^2
FPC_u_MeV         = 931.49432          #unified atomic mass unit (u) (28) MeV/c^2
FPC_u_kg          = 1.6605402e-27      #unified atomic mass unit (u) (10) kg
FPC_EPSILON_0     = 8.854187817e-12    #permittivity of free space F/m
FPC_MU_0          = 1.2566370614e-06   #permeability of free space N/A^2
FPC_ALPHA         = 0.0072973530796448 #fine-structure constant (61)
FPC_r_e           = 2.81794092e-15     #classical electron radius (38) m
FPC_LAMBDA_BAR_e  = 3.86159323e-13     #electron Compton wavelength (35) m
FPC_a_0           = 5.29177249e-11     #Bohr radius (mnucleus= infty) (24) m
FPC_LAMBDA_1EV    = 1.23984244e-06     #wavelength of 1 eV/c particle (37) m
FPC_R_INFINITY_EV = 13.6056981         #Rydberg energy (mnucleus = infinity) (40) eV
FPC_SIGMA_0_BARN  = 0.66524616         #Thomson cross section (18) barn
FPC_MU_B_MeV_T    = 5.78838263e-11     #Bohr magneton (52)  MeV/T
FPC_MU_N_MeV_T    = 3.15245166e-14     #nuclear magneton (28) MeV/T
FPC_E_M_e         = 175881962000       #electron cyclotron freq./field (53) C/kg (rad/sT)
FPC_E_M_P         = 95788309           #proton cyclotron freq./field (29) C/kg (rad/sT)
FPC_G_SI          = 6.67259e-11        #gravitational constant (85) m^3/kgs^2
FPC_G_P           = 6.70711e-39        #gravitational constant (86) h_bar c (GeV/c^2)^{-2}
FPC_g             = 9.80665            #standard grav. accel., sea level m/s^2
FPC_N_A           = 6.0221367e+23      #Avogadro constant (36)  /mole
FPC_K_B           = 1.380658e-23       #Boltzmann constant (12) J/K
FPC_K_B_EV        = 8.617385e-05       #Boltzmann constant (73) eV/K
FPC_V_MOLAR       = 0.0224141          #molar volume, ideal gas at STP (19) m^3/mole
FPC_LAMBDAT       = 0.002897756        #Wien displacement law constant (24) m K
FPC_SIGMA_SB      = 5.67051e-08        #Stefan-Boltzmann constant (19) W/m^2K^4
FPC_G_F           = 1.16639e-05        #Fermi coupling constant (2)  GeV^{-2}
FPC_SIN2_THETA_W  = 0.23192            #weak mixing angle  (5)  at M_Z
FPC_M_W           = 80.22              #W boson mass (26) GeV/c^2
FPC_M_Z0          = 91.187             #Z_0 boson mass (7) GeV/c^2
FPC_G_S           = 0.117              #strong coupling constant (5))  at M_Z
AUCHARGE          = 1.60217733e-19     #1 a.u. charge in SI, elementary charge
AUMASS            = 9.1093897e-31      #1 a.u. mass in SI, electron rest mass
AUENERGY          = 4.3597482e-18      #1 a.u. energy in SI, hartree energy
AULENGTH          = 5.29177249e-11     #1 a.u. length in SI, bohrradius
AUACTION          = 1.05457266e-34     #1 a.u. action in SI, Planck constant/(2*PI)
AUTIME            = 2.4188843341e-17   #1 a.u. time in SI
AUFORCE           = 8.2387295e-08      #1 a.u. force in  SI
AUVELOCITY        = 2187691.42         #1 a.u. velocity in SI
AUMOMENTUM        = 1.9928534e-34      #1 a.u. momentum in SI
AUEFIELD          = 514220820000       #1 a.u. el. field in SI
AUEDIPOLE         = 8.4783579e-30      #1 a.u. el. dipole in SI
AUMFLUX           = 235051.808         #1 a.u. magn. flux in SI
AUMDIPOLE         = 1.85480308e-23     #1 a.u. magn. flux in  SI
ELECTRONVOLT      = 1.60217733e-19     #1 eV charge in SI
ATOMICMASS        = 1.6605402e-27      #1 amu mass in SI
ASTRONOMICALUNIT  = 149597870000       #1 AE length in SI
MOL               = 6.0221367e+23      #Avogadro constant
au2an             = 0.529177249        #[a.u.] -> [Angstroem]
an2au             = 1.8897259885789    #[Angstroem] -> [a.u.]
au2ev             = 27.211396131788    #[a.u.] -> [eV]
ev2au             = 0.036749308824762  #[eV] -> [a.u.]
au2am             = 0.00054857989586762#[a.u.] -> [atomar mass units]
am2au             = 1822.8885300626    #[atomar mass units] -> [a.u.]
au2fs             = 0.024188843341     #[a.u.] -> [fs]
fs2au             = 41.341373206755    #[fs] -> [a.u.]
ev2kc             = 23.05              #[eV] -> [kcal/mol]
kc2ev             = 0.043364476254697  #[kcal/mol] -> [eV]
au2kc             = 627.50431878767    #[a.u.] -> [kcal/mol]
kc2au             = 0.0015936145299079 #[kcal/mol] -> [a.u.]
au2ic             = 219474.7           #[au] -> [1/cm]
ic2au             = 4.5563338279993e-06#[1/cm] -> [au]
ev2ic             = 8065.54353         #[eV] -> [1/cm]
ic2ev             = 0.0001239842047    #[1/cm] -> [eV]
ev2kj             = 96.485             #[eV] -> [kJoule/mol]
kj2ev             = 0.0103643          #[kJoule/mol] -> [eV]
in2cm             = 2.54               #[inch] -> [cm]
ft2m              = 0.3048             #[feet] -> [cm]
yd2m              = 0.9144             #[yard] -> [cm]
cm2in             = 0.39370078740157   #[inch] -> [cm]
m2ft              = 3.2808398950131    #[feet] -> [cm]
m2yd              = 1.0936132983377    #[yard] -> [cm]

au2va             = 51.4220652         #[au of E-Field] -> [V/A]
va2au             = 1.0/au2va          #[V/A] -> [au of E-Field]

def wcm22va(I):
    return 2.744e-7 * I**0.5

def va2wcm2(E):
    return (E/2.744e-7)**2

# Dictionary of element names and atomic number:
ptable = {
        'H'  :  1,
        'He' :  2,
        'Li' :  3,
        'Be' :  4,
        'B'  :  5,
        'C'  :  6,
        'N'  :  7,
        'O'  :  8,
        'F'  :  9,
        'Ne' : 10}

def factor(unit,physical_quantity):
    try:
        if physical_quantity is 'energy' and unit is 'cm-1':
            return ic2au
        if physical_quantity is 'energy' and unit is 'ev':
            return ev2au
        if physical_quantity is 'time' and unit is 'fs':
            return fs2au
        if physical_quantity is 'efield' and unit is 'va':
            return va2au
    except:
        raise Exception('Conversion of '+physical_quantity+' from '+unit+' not available')



