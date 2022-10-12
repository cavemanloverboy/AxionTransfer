import yt
import numpy as np
from nbodykit.lab import *
import matplotlib.pyplot as plt


# ds = yt.load("IC1024_2000/DD0291/DD0291")
ds = yt.load("JedICs_199/RD0003/RD0003")
BOXSIZE = float(ds.length_unit)
Nx = ds.domain_dimensions[0]


# cg = ds.covering_grid(level=L, left_edge=ds.domain_left_edge, dims=[Nx]*3)

p = yt.AxisProjectionPlot(ds, "z",  ("deposit", "nbody_mass"), ("gas", "density"), width=(BOXSIZE, "kpc/h"))
p.save()



def calculate_power(ds):
    data = ds.all_data()
    px =  np.array([data[('all', 'particle_position_x')], data[('all', 'particle_position_y')], data[('all', 'particle_position_z')]]).T

    N = len(px)
    Data = np.empty(N, dtype=[('Position', ('f8', 3)), ('Mass', 'f8')])
    Data['Position'] = px
    Data['Mass'] = np.ones(N)

    catalog = ArrayCatalog(Data)

    catalog['Position'] *= BOXSIZE # re-normalize units of Position


    mesh = catalog.to_mesh(Nmesh=Nx, BoxSize=BOXSIZE)
    #mesh.save('mesh.bigfile')
    #plt.figure()
    #plt.imshow(mesh.preview(axes=[0,1], Nmesh=512))
    #plt.savefig('mesh.png')

    result = FFTPower(mesh, mode='1d')
    #result.save("power-result.json")
    results = result.run()
    k     = results[0]['k']
    power = results[0]['power']
    # modes = results[0]['modes']

    return k, power


def get_initial_power():
    result = np.genfromtxt("tf199.dat", delimiter=" ", dtype=float).T
    return result[:2]
k, power = calculate_power(ds)
ki, pi = get_initial_power()
import pdb; pdb.set_trace()
plt.figure()
plt.loglog(k, power, label="z=0")
plt.loglog(ki, pi, label="z=199")
plt.legend()
plt.savefig(f"{ds.basename}_power.png")