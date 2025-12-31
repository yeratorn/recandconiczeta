import utils as ut
import math
import numpy as np


def HJ(p,n):
  M = [n,p]
  while M[-1]>0:
    M.append((-M[-2])%M[-1])
  return M

def HJI(p,n):
  H = HJ(ut.inv(p,n),n)
  H.reverse()
  return H

def gcd2dcone(V1,V2):
    d = len(V1)
    if len(V2) != d:
        return 'dimension error'
    final = ()
    for i in range(d):
        for j in range(d):
            if i != j:
                det = V1[i]*V2[j]-V1[j]*V2[i]
                final += (det,)
    return math.gcd(*final)

#Give the (p,n) type of the cone
def type2dcone(V1,V2):
    mult = gcd2dcone(V1,V2)
    if mult == 1:
        print('nonsingular')
        return 0,1 
    for i in range(0,mult):
        vec = i*np.array(V1)+np.array(V2)
        if not (vec%mult).any() and i!=0:
            print(i,mult)
            return i,mult

def res2dcone(V1,V2):
    p,n = type2dcone(V1,V2)
    h = HJ(p,n)
    hi = HJI(p,n)
    s = len(h)-1 
    final = []
    for i in range(s+1):
        vec = (1/n)*(h[i]*np.array(V1)+hi[i]*np.array(V2))
        final.append(list(vec))
    return final 