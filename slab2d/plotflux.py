#
#

from boutdata import collect

from numpy import zeros, mean, sqrt, arange, log, exp

import matplotlib.pyplot as plt

from boututils.linear_regression import linear_regression

path="."

xsep = 156

phi = collect("phi",path=path)
pe = collect("pe", path=path)
tnorm = collect("Tnorm")
nnorm = collect("Nnorm")

pe *= 1.602e-19 * nnorm * tnorm # Convert to Pascals

Bxy = 0.5   # Magnetic field [T]

phi = phi * tnorm  # Convert phi to Volts

# Domain size
dx = collect("dx", path=path)[0,0]
dz = collect("dz", path=path)
g_11 = collect("g_11", path=path)[0,0]
g_33 = collect("g_33", path=path)[0,0]
rho_s0   = collect("rho_s0", path=path)


dr = dx * sqrt(g_11)   # Radial grid spacing, in rho_s
dr *= rho_s0 * 1e3 # Convert to mm

dz = dz * sqrt(g_33) * rho_s0 * 1e3  # Binormal grid spacing, in mm

nt,nx,ny,nz = phi.shape # Get array sizes

rarr = arange(nx) * dr
rsep = rarr[156]

# Differentiate phi in Z to get a radial velocity

dphidz = zeros(phi.shape)
for z in range(nz):
    zp = (z + 1) % nz
    zm = (z - 1 + nz) % nz
    dphidz[:,:,:,z] = 1000.*(phi[:,:,:,zp] - phi[:,:,:,zm])/(2.*dz*Bxy)

# dphidz now contains the ExB velocity in m/s

vp = dphidz.copy()
vp[dphidz < 0.0] = 0.0   # Zero vp where velocity is negative

vm = dphidz.copy()
vm[dphidz > 0.0] = 0.0   # Zero vm where velocity is positive

# Calculate positive and negative fluxes
flux_p = 1.5*pe * vp   # Now has units of J/m^2/s
flux_m = 1.5*pe * vm

# Average over time and z
flux_p_av = mean(mean(flux_p[:,:,0,:], axis=0), axis=1)
flux_m_av = mean(mean(flux_m[:,:,0,:], axis=0), axis=1)

net_flux = flux_p_av + flux_m_av

target_flux = zeros(nx)
for x in range(1,nx-1):
    target_flux[x] = (net_flux[x+1] - net_flux[x-1])/(2.*dr*1e-3)

# Fit

xcoord = rarr - rsep
xin = 30
xout = 200
ai, bi = linear_regression(xcoord[(xsep+2):(xsep+xin)], log(-target_flux[(xsep+2):(xsep+xin)]))
ao, bo = linear_regression(xcoord[(xsep+xin):(xsep+xout)], log(-target_flux[(xsep+xin):(xsep+xout)]))

lambda_qin = -1./bi
lambda_qout = -1./bo

fig, ax1 = plt.subplots()

ax1.plot(rarr-rsep, flux_p_av/1000., '-k', label="Outward")
ax1.plot(rarr-rsep, -flux_m_av/1000., '--r', label="Inward")

ax1.plot([0,0],[0,160],'-k', linewidth=2)
ax1.set_ylim([0,160])
ax1.set_xlim([-70,150])

ax1.set_ylabel(r'Radial flux [kW/m$^2$]')
ax1.set_xlabel("Radius from separatrix [mm]")

# Create second y axis
ax2 = ax1.twinx()

ax2.plot((rarr-rsep)[xsep:-10], -target_flux[xsep:-10]/1e6, label="Target heat flux")
ax2.plot(xcoord[xsep:(xsep+xin)], exp(ai+bi*xcoord[xsep:(xsep+xin)])/1e6, '--k',label=r'$\lambda_q^{\mathrm{near}} = %.1f$ mm' % lambda_qin)
ax2.plot(xcoord[(xsep+xin):(xsep+xout)], exp(ao+bo*xcoord[(xsep+xin):(xsep+xout)])/1e6, '-k',label=r'$\lambda_q^{\mathrm{far}} = %.1f$ mm' % lambda_qout)

ax2.set_ylabel(r"Flux to target [MW/m$^3$]")
ax2.set_yscale("log")
ax2.set_xlim([-70,150])

ax1.legend(loc="upper left")
ax2.legend(loc="upper right")
plt.savefig("flux.pdf")

with open("flux.csv", "w") as f:
  f.write("Radius from separatrix [mm], Outward flux [W/m^2], Inward flux [ARB]\n")
  for r, p, m in zip(rarr-rsep, flux_p_av, -flux_m_av):
    f.write(str(r)+", "+ str(p) + ", " + str(m) + "\n")

plt.show()
