#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 15:29:43 2018

@author: alessandrocoretti
"""

import numpy as np
from numpy import linalg as LA
from scipy.stats import ortho_group
 
N = 40
Nclusters = 5
scale = 100
avg = [(i+1)**6*scale for i in range(Nclusters)]
sigma = np.ones(Nclusters)
diagonal = np.array([np.random.normal(size=(N/Nclusters), loc=avg[i], scale=sigma[i]) for i in range(Nclusters)]).flatten()

scale_dd = 100
scale_scn = 10

A = np.random.rand(N,N)
A = (A + A.T)/2
Ainv = LA.inv(A)
wA, vA = LA.eigh(A)
kA = LA.cond(A)

B = np.random.normal(size=(N,N))
B = (B + B.T)/2
Binv = LA.inv(B)
wB, vB = LA.eigh(B)
kB = LA.cond(B)

I = np.identity(N)
Irand = I*np.sort(diagonal)
#Irand = Irand*scale_dd
wrand, vrand = LA.eig(Irand)

C = A + scale_dd*I
C = C + scale_scn*Irand
Cinv = LA.inv(C)
wC, vC = LA.eigh(C)
kC = LA.cond(C)

#randmat = np.random.rand(N,N)
#randmat = rvs(N)
randmat = ortho_group.rvs(N)
#randmat = randmat*(np.ones((N,N)) + 99*np.identity(N))
randmatinv = LA.inv(randmat)
D = np.matmul(randmatinv, np.matmul(Irand, randmat))
Dsymm = .5*(D+D.T)
wD, vD = LA.eig(D)
wD = np.sort(wD)
wDsymm, vDsymm = LA.eig(Dsymm)

b = np.random.rand(N)
x0 = np.random.rand(N)

np.savetxt("A.inpt", C, delimiter='\t')
np.savetxt("b.inpt", b, delimiter='\t')
np.savetxt("x0.inpt", x0, delimiter='\t')

print LA.cond(D)