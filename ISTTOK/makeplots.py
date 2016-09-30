#!/usr/bin/env python

from boututils.datafile import DataFile
import matplotlib.pyplot as plt
from boutdata import collect
from numpy import linspace, concatenate, ravel

with DataFile("isttok-10ev.nc") as d:
    Rxyz = d["Rxyz"]
    Zxyz = d["Zxyz"]

dirs = ["lowres-5ev/part-01", 
        "lowres-5ev/part-02",
        "lowres-10ev/part-01", 
        "lowres-10ev/part-02"]


jxz = collect("jpar", path=dirs[-1], tind=-1, yind=8)
pxz = collect("pe", path=dirs[-1], tind=-1, yind=8)

nnorm = collect("nnorm")
tnorm = collect("tnorm")
cs0 = collect("cs0")
pnorm = 1.602e-19 * nnorm * tnorm  # Pa
jnorm = 1.602e-19 * nnorm * cs0 * 1e-3 # kA/m^2

levels = linspace(-60, 60, 31)
plt.contourf(Rxyz[:,8,:], Zxyz[:,8,:], jxz[0,:,0,:]*jnorm,levels)
plt.xlabel("Major radius [m]")
plt.ylabel("Height [m]")
plt.title(r"Parallel current [kA/m$^2$]")
plt.colorbar()
plt.savefig("isttok-10ev-jpar.pdf")

plt.figure()
plt.contourf(Rxyz[:,8,:], Zxyz[:,8,:], pxz[0,:,0,:]*pnorm,30)
plt.xlabel("Major radius [m]")
plt.ylabel("Height [m]")
plt.title(r"Plasma pressure [Pa]")
plt.colorbar()
plt.savefig("isttok-10ev-pxz.pdf")
plt.show()

def timeseries(var, x,y,z):
    first = True
    for d in dirs:
        tmp = collect(var, xind=x, yind=y, zind=z, path=d)
        tmp = ravel(tmp)
        if first:
            data = tmp
            first = False
        else:
            data = concatenate([data, tmp])
    return data
    

time = timeseries("t_array", 0,0,0)
jout = timeseries("jpar", 10, 8, 26) # y = 8, z=26 : outboard midplane
jin = timeseries("jpar", 10, 8, 154) # y = 8, z = 154 : Inboard midplane

wci = collect("Omega_ci")
time *= 1e3/wci # ms

plt.figure()
plt.plot(time, jout*jnorm, label="Outboard midplane")
plt.plot(time, jin*jnorm, label="Inboard midplane")
plt.xlabel("Time [ms]")
plt.legend()
plt.ylabel(r"$J_{||}$ [kA/m$^2$]")
plt.savefig("isttok-10ev-jpar-time.pdf")

#############################

pout = timeseries("Pe", 10, 8, 26) # y = 8, z=26 : outboard midplane
pin = timeseries("Pe", 10, 8, 154) # y = 8, z = 154 : Inboard midplane

plt.figure()
plt.plot(time, pout, label="Outboard midplane")
plt.plot(time, pin, label="Inboard midplane")
plt.xlabel("Time [ms]")
plt.legend()
plt.ylabel(r"Pressure [Pa]")
plt.savefig("isttok-10ev-pe-time.pdf")

######################## SOL fluctuations

pout = timeseries("Pe", 30, 8, 26) # y = 8, z=26 : outboard midplane
pin = timeseries("Pe", 30, 8, 154) # y = 8, z = 154 : Inboard midplane

plt.figure()
plt.plot(time, pout, label="Outboard midplane")
plt.plot(time, pin, label="Inboard midplane")
plt.xlabel("Time [ms]")
plt.legend()
plt.ylabel(r"Pressure [Pa]")
plt.savefig("isttok-10ev-pe-time-sol.pdf")
plt.show()
