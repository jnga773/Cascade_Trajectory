#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 14:05:49 2018

@author: jnga773
"""

import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rcParams
rcParams['text.usetex'] = True
rcParams['text.latex.unicode'] = True
rcParams['font.family'] = 'serif'

plt.close("all")

params = ""
direct = "./data_files/" + params

# Extract system variables
omega = float(open(direct + "tau.txt").readlines()[0])
delta = float(open(direct + "tau.txt").readlines()[1])
xi    = float(open(direct + "tau.txt").readlines()[2])
alpha = float(open(direct + "tau.txt").readlines()[3])
D_a   = float(open(direct + "tau.txt").readlines()[4])
kappa_a = float(open(direct + "tau.txt").readlines()[5])
D_b   = float(open(direct + "tau.txt").readlines()[6])
kappa_b = float(open(direct + "tau.txt").readlines()[7])
gamma = 1.0

# time values
tau = np.loadtxt(fname=direct + "tau.txt", dtype='float', usecols=(0,), skiprows=9)
dt = tau[1] - tau[0]

# Extract correlation and cross correlation value
corr = np.genfromtxt(fname=direct + "correlation.txt", dtype='float', usecols=(0,))
cross_corr = np.genfromtxt(fname=direct + "correlation.txt", dtype='float', usecols=(1,))

fig, ax = plt.subplots(figsize=[8,5])

left, bottom, width, height = [0.5, 0.3, 0.4, 0.4]
#ax2 = fig.add_axes([left, bottom, width, height])

ax.plot(tau, corr, label=r'Correlation from a, $\Delta_{a}=%s$'%(D_a))
ax.plot(tau, cross_corr, label=r'Cross correlation from b, $\Delta_{b}=%s$'%(D_b))

ax.set_xlabel(r'$\gamma \tau$')
ax.set_ylabel(r'$g^{(2)}(\tau)$')
ax.legend()

#ax2.plot(tau, corr)
#ax2.set_xlim(-0.05, 3)
#ax2.set_xlabel(r'$\gamma\tau$')
#ax2.set_ylabel(r'$g^{(2)}(\tau)$')
#
#ax.set_title(r'$\Omega=%s, \delta=%s, \xi=%s, \alpha=%s, \kappa_{a}=%s, \kappa_{b}=%s$'%(omega, delta, xi, alpha, kappa_a, kappa_b))
fig.tight_layout()

#fig.savefig("../Images/traj/omega=%s_D_a=%s_D_b=%s_kappa_a=%s_kappa_b=%s.pdf"%(omega, D_a, D_b, kappa_a, kappa_b))