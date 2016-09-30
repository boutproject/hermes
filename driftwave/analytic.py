#
#

from numpy import sqrt

def dispersion(kperprho, sparn, tau):
    """
    Input
    
    kperprho = k_\perp * rho_s
    sparn    = sigma_|| / w_star
    tau      = T_i / T_e
    
    Returns
    
    Complex frequency / omega_star
    
    """
    
    a = 1.
    b = (tau + (1. + tau)*1j*kperprho**2*sparn) + 1j*sparn
    c = -1j*sparn
    
    w1 = (-b + sqrt(b**2 - 4.*a*c) ) / (2.*a)
    w2 = (-b - sqrt(b**2 - 4.*a*c) ) / (2.*a)
    
    return w1, w2


if __name__ == "__main__":
    from numpy import linspace
    import matplotlib.pyplot as plt
    
    kpar = 10.**linspace(-3,1., 100)
    tau = 1.
    kperp = 1.
    #for kperp in [1e-4, 0.1, 1.0, 10.]:
    for tau,col in [(0.0, 'r'), 
                    (1.0, 'g'),
                    (2.0, 'b')]:
        w1, w2 = dispersion(kperp, kpar, tau)
        
        plt.plot(kpar, w1.imag, label="Growth-rate tau="+str(tau), color=col)
        plt.plot(kpar, w1.real, label="Frequency tau="+str(tau), linestyle='--', color=col)
        
    plt.xscale('log')
    plt.legend(loc='upper left')
    plt.grid()
    plt.show()
    
