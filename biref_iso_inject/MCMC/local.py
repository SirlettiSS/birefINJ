
# general
import numpy as np


# physical units
Tcmb  = 2.726e6       # CMB temperature
ac2rad = np.pi/10800. # arcmin -> rad
eV2Mpc = 806554.815354305*3.0856775814914e22*(2*np.pi)

# Planck or CamSpec NPIPE + CMB Lensinging + DESI DR2 BAO parameters
Planck = False


if Planck:
    # cosmological parameters (Planck 2018, arXiv:1807.06209, TT+TE+EE+lowE+lensing)
    H0    = 67.36
    Om    = 0.3153
    As    = 2.1e-09
    ns    = 0.9649
    ombh2 = 0.02237
    zcmb  = 1059.94
    zeq   = 3402

else:
    # cosmological parameters (Toshiya's best-fit: CamSpec NPIPE + CMB Lensing + DESI DR2 BAO)
    H0    = 68.465
    Om    = 0.29917
    As    = 2.25614441e-09
    ns    = 0.97006
    ombh2 = 0.02239
    zcmb  = 1089.67 
    zeq   = 3402    

# derived cosmological parameters
acmb  = 1./(1.+zcmb)
aeq   = 1./(1.+zeq)
h0    = H0/100.
#H0Mpc = H0/c.C
H0Mpc = H0/299792.458
Ob    = ombh2/h0**2
Or    = Om*aeq
Ov    = 1-Om # flat universe
omch2 = Om*h0**2 - ombh2
cps   = {'H0':H0,'Om':Om,'Ov':Ov,'w0':-1,'wa':0.}


# phi parameters
def params_phi(model='ALP'):
    '''
    n      : n=1 for ALP and n=2 for EDE potential
    tf     : f_a/phi_ini; the symmetry breaking scale normalized by phi_ini, dimensionless
    '''
    if model=='ALP':
        n, tf = 1, 100.
    if model=='EDE':
        n, tf = 2, np.sqrt(8*np.pi)
    return n, tf


def Ea(a,rad=1.): # H(a)/H0
    return np.sqrt( Om/a**3 + Ov + rad*Or/a**4 )


def dlnEada(a):
    return -(1.5*a*Om+2*Or)/(a**5*Ea(a)**2)



def V_phi(x,logm,tf,n,deriv=0):
    '''
    x        : pseudoscalar field normalized by phi_ini, dimensionless
    logm     : logarithm of the mass parameter in the unit of eV
    tf       : the symmetry breaking scale normalized by phi_ini, dimensionless
    deriv    : the order of derivative with respect to x, i.e. the normalized field phi/phi_ini
               (not with respect to the dimensional field phi)
    
    V_phi(x) : renormalized potential wrt phi_ini
    '''
    
    m = 10**(logm)*eV2Mpc
    
    if deriv == 0: # potential V_phi(x)
        return m**2*tf**2*(1.-np.cos(x/tf))**n

    if deriv == 1: # dV_phi/dx, equal to (1/phi_ini) dV/dphi since V_phi is already normalized wrt phi_ini
        return n*m**2*tf*(1.-np.cos(x/tf))**(n-1)*np.sin(x/tf)

    if deriv == 2: # d^2V_phi/dx^2, equal to (1/phi_ini^2) d^2V/dphi^2
        if n == 1:
            return m**2*np.cos(x/tf)
        else:
            return n*m**2*(1.-np.cos(x/tf))**(n-2)*(-1.+np.cos(x/tf)+n*np.sin(x/tf)**2)


def EoM_phi(a, x, logm, tf, n, xi=None, eps=1e-3, z_start=1e2, z_end=10, phi_ini=1.0):
    
    '''
    # Solving the Klein-Gordon equation in presence of an energy injection term,
    # written as a first-order system in the expansion variable a (not conformal time).
    
    #   dx0/da = x1
    #   dx1/da = d^2x0/da^2 = - ( 4/a + dlnH/da ) dx0/da - V_phi(x0, deriv=1)/calH**2 + Src/calH**2
    # where x0 is varphi = phi/phi_ini and x1 is dvarphi/da.
    # V_phi(x0, deriv=1) is already normalized wrt phi_ini, since it is defined
    # in terms of x = phi/phi_ini and tf = f/phi_ini, so V_phi(deriv=1) = (1/phi_ini) dV/dphi.
    
    # Src is the regularized energy injection term:
    #   Src = Q / (dvarphi/da)   with   Q = xi * calH * rho_M  (model for the injected energy density rate)
    # regularized as Q * x1 / (x1**2 + eps**2) to avoid the singularity at dvarphi/da = 0
    # (field turning points), while preserving the correct sign via the x1 numerator.
    '''
    
    calH = a * H0Mpc * Ea(a) # conformal expansion rate
    m = 10**(logm) * eV2Mpc  # mass
    z = 1./a - 1.
    
    if xi is not None and z > z_end and z < z_start:
        frac = 1e-3  # rho_m0/rho_phi_ini, i.e. energy density fraction of ALP to matter
        V_ini = V_phi(phi_ini, logm, tf, n, deriv=0)
        Q    = (xi/a) * (V_ini / frac / a**3) # energy injection rate to the ALP field
        Src  = Q * x[1] / (x[1]**2 + eps**2)  # regularized Q/(dvarphi/da)
        
        
    else:
        Src = 0.
    
    return [ x[1], - ( (4./a + dlnEada(a))*x[1] + ( V_phi(x[0], logm, tf, n, deriv=1) - Src ) / calH**2 ) ]


