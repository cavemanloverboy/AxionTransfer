
import sys, platform, os
import matplotlib
from matplotlib import pyplot as plt
import numpy as np

import camb
from camb import model, initialpower

#Set up a new set of parameters for CAMB
pars = camb.CAMBparams()

# Planck 2018
Omega_m         = 0.3111  #checked
Omega_L         = 0.6889  #checked
Omega_b         = 0.049   #checked
H0              = 67.66   #checked
sigma_8         = 0.8102  #checked
nspec           = 0.9665  #checked

pars.set_cosmology(
    H0=H0, 
    ombh2=Omega_b*(H0/100)**2, 
    omch2=(Omega_m-Omega_b)*(H0/100)**2, 
    mnu=0.06, 
    omk=0, 
    tau=0.0523,
)
pars.InitPower.set_params(
    As=2.09324e-9, # Amplitude
    ns=nspec, # Spectral Index
    r=0, # Tensor to scalar
)
results = camb.get_results(pars)

#Linear spectra
z0 = 9999
redshifts = [z0, 3]
redshifts.sort()
pars.set_matter_power(redshifts=redshifts, kmax=1000.0)
pars.NonLinear = model.NonLinear_none
results = camb.get_results(pars)
s8 = np.array(results.get_sigma8())


plt.figure()
kh, z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=681.52510591, npoints = 1000)
plt.loglog(kh, pk[-1], label=str(z[-1]))
plt.xlabel("k [h / Mpc]")
plt.ylabel("P(k)")
plt.grid()
plt.savefig(f"power{z0}.png")

kh, z, pk = results.get_matter_power_spectrum(minkh=12.11942063, maxkh=681.52510591, npoints = 1000)

fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True,  gridspec_kw={'height_ratios': [3, 1]})
for i in range(len(pk)):
    axes[0].loglog(kh, pk[i], label=str(int(z[i])))


# Jed Transfer
jed_ktk = np.genfromtxt("jed_transfer.csv", delimiter=",")

# Unpack log10_ktilde and transfer
jed_log10_ktilde, jed_tk = jed_ktk.T

# Convert ktilde to k in units of h/Mpc
jed_k = 10**jed_log10_ktilde * 82.0 / (H0 / 100)

# interpolate to common grid
jed_interp_pk = pk[-1] * np.interp(kh, jed_k, jed_tk)

print(kh)
print(jed_k)
print(jed_tk)

axes[0].loglog(kh, jed_interp_pk, label=f"{z0} jed")
axes[0].set_ylabel("P(k)")
axes[1].loglog(kh, np.interp(kh, jed_k, jed_tk))
axes[1].set_ylabel("P(k) Ratio")
axes[0].legend()
axes[0].grid()
axes[1].grid()
axes[1].set_xlabel("k [h / Mpc]")
fig.savefig(f"powers{z0}.png")

# # FROM MUSIC
# ss >> k;
# ss >> Tkc;   // cdm
# ss >> Tkb;   // baryon
# ss >> dummy; // photon
# ss >> dummy; // nu
# ss >> dummy; // mass_nu
# ss >> Tktot; // total
# ss >> dummy; // no_nu
# ss >> dummy; // total_de
# ss >> dummy; // Weyl
# ss >> Tkvc;  //>[150609SH: add] // v_cdm
# ss >> Tkvb;  //>[150609SH: add] // v_b
# ss >> dummy; // v_b-v_cdm
k = kh
Tkc = jed_interp_pk
Tkb = np.zeros_like(Tkc)
dummy = np.zeros_like(Tkc)
dummy2 = dummy
dummy3 = dummy
Tktot = Tkc
dummy4 = dummy
dummy5 = dummy
dummy6 = dummy
tkvc = dummy # new dummy, replacing with 2lpt
tkvb = dummy # new dummy, replacing with 2lpt
dummy7 = dummy


DataOut = np.column_stack((kh, Tkc, Tkb, dummy, dummy2, dummy3, Tktot, dummy4, dummy5, dummy6, tkvc, tkvb, dummy7))
np.savetxt(f'tf{z0}.dat', DataOut)

# # Binary Format
# import struct
# with open('tf.dat', 'w') as your_dat_file:  
#     your_dat_file.write(struct.pack('d'*len(kh), *kh))
#     your_dat_file.write(struct.pack('d'*len(jed_interp_pk), *jed_interp_pk))
#     your_dat_file.write(struct.pack('d'*len(jed_interp_pk), *jed_interp_pk*0))
#     your_dat_file.write(struct.pack('d'*len(jed_interp_pk), *jed_interp_pk*0))
#     your_dat_file.write(struct.pack('d'*len(jed_interp_pk), *jed_interp_pk*0))
#     your_dat_file.write(struct.pack('d'*len(jed_interp_pk), *jed_interp_pk*0))
#     your_dat_file.write(struct.pack('d'*len(jed_interp_pk), *jed_interp_pk*0))
#     your_dat_file.write(struct.pack('d'*len(jed_interp_pk), *jed_interp_pk*0))
#     your_dat_file.write(struct.pack('d'*len(jed_interp_pk), *jed_interp_pk*0))
#     your_dat_file.write(struct.pack('d'*len(jed_interp_pk), *jed_interp_pk*0))
#     your_dat_file.write(struct.pack('d'*len(jed_interp_pk), *jed_interp_pk*0))
#     your_dat_file.write(struct.pack('d'*len(jed_interp_pk), *jed_interp_pk*0))
#     your_dat_file.write(struct.pack('d'*len(jed_interp_pk), *jed_interp_pk*0))

# with open('your_data.txt', 'w') as your_dat_file:
#     for i in range(len(kh)):
#         your_dat_file.write(f"{kh[i]} {jed_interp_pk[i]} 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n")