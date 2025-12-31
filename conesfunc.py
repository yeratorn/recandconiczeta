import numpy as np
import math
import itertools as iter

maxN = 100

def quot(n,q):
    if len(n) != len(q):
        return 'dimension error'
    prod = 1
    for i in range(len(n)):
        if n[i] == 0 and q[i] == 0:
            prod = prod
        elif n[i] == 0 and q[i] != 0:
            prod = 0
        else:
            prod *= n[i]**(-q[i])
    return prod

def cone(L,N = maxN):
  dim = len(L[0])
  for l in L:
    if len(l) != dim:
      return 'dimension error'
  r = len(L)
  lis = list(iter.product(range(1,N),repeat = r))
  final = []
  for N in lis:
    vec = np.zeros(dim)
    for i in range(r):
      vec += N[i]*np.array(L[i])
    final.append(vec.astype(float).tolist())
  return final


def zetacone(L,q):
    C = cone(L)
    final = 0
    for P in C:
        #quot(P,q)
        final += quot(P,q)
    return final


def conegcd(L,N = maxN):
  dim = len(L[0])
  for l in L:
    if len(l) != dim:
      return 'dimension error'
  r = len(L)
  lis = iter.product(range(1,N),repeat = r)
  final = []
  for cf in lis:
      if math.gcd(*cf)==1:
        vec = np.zeros(dim)
        for i in range(r):
          vec += cf[i]*np.array(L[i])
          newvec = vec.astype(float).tolist()
        final.append(newvec)
  return final   

def zetaconegcd(L,q):
    C = conegcd(L)
    final = 0
    for P in C:
        final += quot(P,q)
    return final