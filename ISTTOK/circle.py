# Generates an input mesh for circular, large aspect-ratio
# simulations:
#
# o Constant magnetic field
# o Curvature output as a 3D logB variable
# o Z is poloidal direction
# o Y is parallel (toroidal)
#
# NOTE: This reverses the standard BOUT/BOUT++ convention
#       so here Bt and Bp are reversed
#

from __future__ import division

from numpy import zeros, ones, ndarray, pi, cos, sin, outer, linspace,sqrt,exp, arange

from boututils.datafile import DataFile # Wrapper around NetCDF4 libraries

def generate(nx, ny, nz,
             R = 2.0, r=0.2, # Major & minor radius
             dr=0.05,  # Radial width of domain
             Bt=1.0,   # Toroidal magnetic field
             q=5.0,    # Safety factor
             mxg=2,    # Guard cells in X
             file="circle.nc",   # Output file name
             ixseps=None,    #
             ):
    """
    Inputs
    ------
    nx, ny, nz  - Size of the domain, including X guard cells
    
    R       - Major radius [m]
    r       - minor radius [m]
    
    Bt      - Toroidal magnetic field [T]
    q       - Safety factor  (constant)
    
    mxg     - Number of points in X guard cell
    
    file    - Output file name
    
    ixseps  - Separatrix location; number of points in core
    
    
    """
    
    if ixseps == None:
        ixseps = nx   # Default is all points in the core
    
    # q = rBt / RBp
    Bp = r*Bt / (R*q)
    B = sqrt(Bt**2 + Bp**2)
    
    # Minor radius as function of x. Choose so boundary
    # is half-way between grid points

    h = dr / (nx - 2.*mxg) # Grid spacing in r
    rminor = linspace(r - 0.5*dr - (mxg-0.5)*h,
                      r + 0.5*dr + (mxg-0.5)*h,
                      nx)
    
    # mesh spacing in x and y
    dx = ndarray([nx,ny])
    dx[:,:] = r*Bt*h      # NOTE: dx is toroidal flux

    dy = ndarray([nx,ny])
    dy[:,:] = 2.*pi / ny

    # Calculate poloidal angle, R and Z for plotting
    theta = ndarray([nx,ny,nz])
    Rxyz = ndarray([nx, ny, nz])
    Zxyz = ndarray([nx, ny, nz])

    logB = ndarray([nx, ny, nz])
    
    for y in range(ny):
        dtheta = y * 2.*pi / ny / q # Change in poloidal angle
        # Note: Z is -theta
        thetaz = linspace(0,-2.*pi, nz, endpoint=False) + dtheta
        
        theta[:,y,:] = outer( ones(nx), thetaz )
        Rxyz[:,y,:] = R + outer(rminor, cos(thetaz))
        Zxyz[:,y,:] = outer(rminor, sin(thetaz))
                              
    
        # LogB = log(1/(1+r/R cos(theta))) =(approx) -(r/R)*cos(theta)
        logB = -(1./R)*outer(rminor, cos(thetaz))
    
    # Curvature as 3D vector components
    bxcvx = -2.*Bt**2 * r / (B**2 * R) * sin(theta)
    bxcvz = 2.*Bt/(B**2 * r * R) * cos(theta)
    
    # Shift angle from one end of y to the other
    ShiftAngle = ndarray([nx])
    ShiftAngle[:] = -2.*pi / q
    
    Rxy = ndarray([nx,ny])
    Rxy[:,:] = r    # NOTE  : opposite to standard BOUT convention

    Btxy = ndarray([nx,ny])
    Btxy[:,:] = Bp

    Bpxy = ndarray([nx,ny])
    Bpxy[:,:] = Bt

    Bxy = ndarray([nx,ny])
    Bxy[:,:] = sqrt(Bt**2 + Bp**2)

    hthe = ndarray([nx,ny])
    hthe[:,:] = R

    print("Writing to file '"+file+"'")

    f = DataFile()
    f.open(file, create=True)

    # Mesh size
    f.write("nx", nx)
    f.write("ny", ny)
    f.write("nz", nz)

    # Mesh spacing
    f.write("dx", dx)
    f.write("dy", dy)

    # Metric components
    f.write("Rxy", Rxy)
    f.write("Btxy", Btxy)
    f.write("Bpxy", Bpxy)
    f.write("Bxy", Bxy)
    f.write("hthe", hthe)

    f.write("ixseps1", ixseps)

    # Shift
    f.write("ShiftAngle", ShiftAngle);

    # Curvature
    f.write("logB", logB)
    f.write("bxcvx", bxcvx)
    f.write("bxcvz", bxcvz)
    
    
    # Input parameters
    f.write("R", R)
    f.write("r", r)
    f.write("dr", dr)
    f.write("Bt", Bt)
    f.write("q", q)
    f.write("mxg", mxg)

    f.write("theta", theta)
    f.write("Rxyz", Rxyz)
    f.write("Zxyz", Zxyz)

    f.close()

# Modified Tanh profile
# Formula from R.J.Groebner Nucl. Fusion 41, 1789 (2001)

def mtanh(alpha, z):
    return ( (1 + alpha*z)*exp(z) - exp(-z) ) / ( exp(z) + exp(-z) )

def mtanh_profile(x, xsym, hwid, offset, pedestal, alpha):
    """
    Profile with linear slope in the core
    
    xsym - Symmetry location
    hwid - Width of the tanh drop
    offset - zero offset
    pedestal - Height of the tanh 
    """
    b = 0.5*(offset + pedestal)
    a = pedestal - b
    z = (xsym - x) / hwid
    return a * mtanh(alpha, z) + b

def profile(filename, name, offset, pedestal, hwid=0.1, alpha=0.1):
    """
    Calculate a radial profile, and add to file
    """
    with DataFile(filename, write=True) as d:
        nx = d["nx"]
        ny = d["ny"]
        x = arange(nx)
        ix = d["ixseps1"]
        
        prof = mtanh_profile(x, ix, hwid*nx, offset, pedestal, alpha)
        
        prof2d = zeros([nx,ny])
        for y in range(ny):
            prof2d[:,y] = prof
            
        #d[name] = prof2d
	d.write(name, prof2d)

def coordinates(nx, ny, nz,
                R = 2.0, r=0.2, # Major & minor radius
                dr=0.05,  # Radial width of domain
                Bt=1.0,   # Toroidal magnetic field
                q=5.0,    # Safety factor
                mxg=2
                ):
    """
    Returns coordinates (R,Z) as a pair of arrays
    
    """

    h = dr / (nx - 2.*mxg) # Grid spacing in r
    rminor = linspace(r - 0.5*dr - (mxg-0.5)*h,
                      r + 0.5*dr + (mxg-0.5)*h,
                      nx)
    
    print("Grid spacing: Lx = %e, Lz = %e" % (h, 2.*pi*r/nz))
    
    Rxyz = ndarray([nx, ny, nz])
    Zxyz = ndarray([nx, ny, nz])
    
    for y in range(0,ny):
        dtheta = y * 2.*pi / ny / q # Change in poloidal angle
        theta = linspace(0,2.*pi, nz, endpoint=False) + dtheta
        
        Rxyz[:,y,:] = R + outer(rminor, cos(theta))
        Zxyz[:,y,:] = outer(rminor, sin(theta))

    return Rxyz, Zxyz



if __name__ == "__main__":
    # Generate ISTTOK mesh

    import os

    # q = 5 test case
    generate(68,16, 256, R=0.46, r=0.085, dr=0.02, Bt=0.5, q=5.0, ixseps=20, file="isttok-68x16x256-q5-5ev.nc")
    profile("isttok-68x16x256-q5-5ev.nc", "Te0", 0.2, 5)
    profile("isttok-68x16x256-q5-5ev.nc", "Ni0", 1e-4, 1e-2)

    # q = 5 high resolution
    generate(132,32, 512, R=0.46, r=0.085, dr=0.02, Bt=0.5, q=5.0, ixseps=20, file="isttok-132x32x512-q5-5ev.nc")
    profile("isttok-132x32x512-q5-5ev.nc", "Te0", 0.2, 5)
    profile("isttok-132x32x512-q5-5ev.nc", "Ni0", 1e-4, 1e-2)
   
    # q = 7 low resolution
    generate(68,16, 256, R=0.46, r=0.085, dr=0.02, Bt=0.5, q=7.0, ixseps=20, file="isttok-68x16x256-q7-5ev.nc")
    profile("isttok-68x16x256-q7-5ev.nc", "Te0", 0.2, 5)
    profile("isttok-68x16x256-q7-5ev.nc", "Ni0", 1e-4, 1e-2)
    
    # Increase temperature to 10eV
    os.system("cp isttok-68x16x256-q7-5ev.nc isttok-68x16x256-q7-10ev.nc")
    profile("isttok-68x16x256-q7-10ev.nc", "Te0", 0.4, 10.)
    
    # 20eV
    os.system("cp isttok-68x16x256-q7-10ev.nc isttok-68x16x256-q7-20ev.nc")
    profile("isttok-68x16x256-q7-20ev.nc", "Te0", 0.8, 20.)

