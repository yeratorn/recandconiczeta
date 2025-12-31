import math
import random
import itertools as iter
from typing import List

def inv(p,n):
  for i in range(1,n):
    if (p*i)%n == 1:
      return i
  
def dot(L1,L2):
    if len(L1) != len(L2):
        return 'dimension error'
    final = 0
    for i in range(len(L1)):
        final += L1[i]*L2[i]
    return final

def hadprod(L1,L2):
    if len(L1) != len(L2):
        return 'dimension error'
    final = []
    for i in range(0,len(L1)):
        final.append(L1[i]*L2[i])
    return final

def zeta(q,N = 1000):
  s = 0
  for n in range(1,N):
      s += n**(-q)
  return s

def U(v,k):
    r = len(v)-1
    lis = list(iter.product([-1,1],repeat = r))
    vv = [v[i] for i in range(len(v)) if i != k]
    U = []
    for u in lis:
        if dot(u,vv)%v[k] == 0: 
            U.append(u)
    return U

def prodlist(L):
    prod = 1
    for i in range(len(L)):
        prod *= L[i]
    return prod

def frac(x):
    return x-math.floor(x)

def Dedekind(v,k,functions):
    r = len(v)-1
    idx = [i for i in range(r+1) if i != k]
    vk = v[k]
    final = 0
    for i in range(1,vk):
        prod = 1
        for j in idx:
            x = frac(v[j]*i/vk)
            prod *= functions[j](x)
        final += prod
    return final

def randvs(r, low=1, high=15, max_tries=10000, rng=None):
    """
    Genera una lista v = [v0,...,vr] de r+1 enteros > 1, coprimos a pares.
    
    Parámetros:
      r: genera r+1 números
      low, high: rango [low, high] para muestrear
      max_tries: tope de intentos antes de fallar
      rng: opcional, un random.Random(...) para reproducibilidad
    
    Retorna:
      lista de longitud r+1 con enteros > 1 y gcd(vi,vj)=1 para i!=j
    """
    if r < 0:
        raise ValueError("r debe ser >= 0")
    if low < 2:
        low = 2
    if high < low:
        raise ValueError("high debe ser >= low")

    rng = rng or random
    v = []
    tries = 0

    while len(v) < r + 1:
        if tries >= max_tries:
            raise RuntimeError(
                f"No se pudo generar {r+1} números coprimos a pares "
                f"en [{low},{high}] tras {max_tries} intentos. "
                "Prueba aumentando 'high' o 'max_tries'."
            )

        x = rng.randint(low, high)
        if all(math.gcd(x, y) == 1 for y in v):
            v.append(x)
        tries += 1

    return v

def cotas(k,v,q):
    vl = 1000
    ind = -1
    prod = 1
    for j in range(0,len(v)):
        if j != k and v[j]<vl:
            vl = v[j]
            ind = j
    for j in range(len(v)):
       if j != k:
           prod *= zeta(q[j])
    #print(vl)
    cotaY =  2*vl**(q[ind])*prod
    cotaZ = prod
    return cotaY, cotaZ, vl**(q[ind])

def S(q: List[int]) -> List[List[int]]:
    """
    Dado q = [q0, ..., qr] (enteros >= 0), retorna la lista de todos los
    s = [s0, ..., sr] tales que 0 <= sj <= qj para todo j.
    """
    if any((not isinstance(x, int)) or x < 0 for x in q):
        raise ValueError("q debe ser una lista de enteros >= 0")

    return [list(s) for s in iter.product(*(range(qj + 1) for qj in q))]