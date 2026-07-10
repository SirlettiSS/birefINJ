def EoM_phi_mass_potential(t, y, mass):
    '''
    EoM equation of mass potential model
    \ddot{\phi} + 3H(t)\dot{\phi}+ m^2\phi = 0

    Args:
    t (float) : proper time 
    y [float, float] : [0th derivative , 1st derivative]
    mass (float) : the mass of ALP (eV)

    Returns:
    [float, float] :  [1st derivative, 2nd derivative] 
    '''

    return [ y[1], - 3 * time2Hubble(t) * y[1] - (eV2MpcInv**2) * (mass**2) * y[0]] 

def EoM_phi_cos_potential(t, y, mass, f):
    '''
    EoM equation of Axion (cos) potential model
    \ddot{\phi} + 3H(t)\dot{\phi}+ m^2*f*\sin(\phi/f)= 0

    Args:
    t (float) : proper time 
    y [float, float] : [0th derivative , 1st derivative]
    mass (float) : the mass of ALP (eV)
    f (float) : the constant of the potential

    Returns:
    [float, float] :  [1st derivative, 2nd derivative] 
    '''
        
    return [ y[1], - 3 * time2Hubble(t) * y[1] - (eV2MpcInv**2) * (mass**2) * f * np.sin(f * y[0])] 
