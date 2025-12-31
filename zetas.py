import numpy as np
from sympy import bernoulli
import utils as ut
import math
import itertools as iter


maxN = 100

def vperp(v,k,N = maxN):
    r = len(v)-1
    C = list(range(-N,N+1))
    C.remove(0)
    lis = iter.product(C,repeat = r)
    #creo generadores
    G = []
    for j in range(r+1):
        g = [0]*(r+1)
        if j != k:
            g[j] = v[k]
            g[k] = -v[j]
            G.append(g)
    final = []
    for n in lis:
        vec = np.zeros(r+1)
        for j in range(0,r):
            ra =  (n[j]/v[k])
            vec += ra*np.array(G[j])
            vec = vec
        val = vec[k]
        if abs(val - round(val)) < 1e-9:
            vec[k] = round(val)     # fuerzo a entero exacto
            final.append(vec)
        #if vec[k] == int(vec[k]):
         #   final.append(vec)
    return final

def bernrec(v,q):
    fs = [lambda x, s=qq: bernoulli(s, x) for qq in q]
    r = len(v)-1
    final = 0
    for k in range(r+1):
        if q[k] == 1:
            final += ut.Dedekind(v, k, fs)
    return final

def zetakvq(k,v,q):
    Vp = vperp(v,k)
    final = 0
    for n in Vp:
        prod = 1
        for j in range(len(n)):
            if j != k:
                prod *= n[j]**(-q[j])
        final += prod
    return final     

def Iqk(q,k):
    Iqk = 1
    r = len(q)-1
    for j in range(0,r+1):
        if j != k:
            Iqk *= ((2j*math.pi)**(-q[j]))*math.factorial(q[j])
    return Iqk

def bernreczeta(v,q):
    r = len(v)-1
    final = 0
    for k in range(0,r+1):
        if q[k] == 1:
            final += Iqk(q,k)*v[k]*zetakvq(k,v,q)  
    return final

def deltabernrec(q):
    power = 0
    factor = 1
    for k in range(len(q)):
        if q[k] == 1:
            power += 1
    if power != 1:
        factor = (1-(-1)**power)*(2**(-power))
    for k in range(len(q)):
        if q[k] != 1:
            factor *= bernoulli(q[k],0)
    return factor

##ZETA VALUES Y and Z

def Zvq(v,q):
    k=0
    Points = vperp(v,k)
    #print(Points)
    final = 0
    for P in Points:
        #print(P)
        quotient = 1
        for j in range(len(v)):
            if P[j] != 0:
                quotient *= P[j]**(-q[j])
            elif P[j] == 0:
                quotient = 0
        final += quotient
    return final

def Zperp(k,v,q):
    Points = vperp(v,k)
    final = 0
    #q[k] = 0
    for P in Points:
        #print(P)
        quotient = 1
        if P[k] == 0:
            for j in range(0,len(v)):
                if j != k:
                    quotient *= P[j]**(-q[j])
            #print(P)
            #final += quot(P,q)
            final += quotient
    return final

def Yperp(k,v,q):
    Points = vperp(v,k)
    #print(Points)
    final = 0
    #q[k] = 0
    for P in Points:
        #print(P)
        quotient = 1
        if P[k] != 0:
            for j in range(0,len(v)):
                if j != k:
                    quotient *= P[j]**(-q[j])
            #print(P)
            #final += quot(P,q)
            final += quotient
    return final

##ORTHANT ZETAS
def vperpu(v,k,u,N = maxN):
    r = len(v)-1
    C = list(range(1,N+1))
    lis = iter.product(C,repeat = r)
    #creo generadores
    G = []
    v = ut.hadprod(v, u)
    for j in range(r+1):
        g = [0]*(r+1)
        if j != k:
            g[j] = v[k]
            g[k] = -v[j]
            G.append(g)
    final = []
    for n in lis:
        vec = np.zeros(r+1)
        for j in range(0,r):
            ra =  (n[j]/v[k])
            vec += ra*np.array(G[j])
            vec = vec
        val = vec[k]
        if abs(val - round(val)) < 1e-9 and round(val)>=0:
            vec[k] = round(val)     # fuerzo a entero exacto
            final.append(vec)
        #if vec[k] == int(vec[k]):
         #   final.append(vec)
    return final

def Yperpu(v,k,q,u):
    Points = vperpu(v,k,u)
    final = 0
    for P in Points:
        #print(P)
        quotient = 1
        if P[k] != 0:
            for j in range(0,3):
                if j != k:
                    quotient *= P[j]**(-q[j])
            #print(P)
            #final += quot(P,q)
            final += quotient
    return final

def Zperpu(v,k,q,u):
    Points = vperpu(v,k,u)
    final = 0
    for P in Points:
        #print(P)
        quotient = 1
        if P[k] == 0:
            for j in range(0,3):
                if j != k:
                    quotient *= P[j]**(-q[j])
            #print(P)
            #final += quot(P,q)
            final += quotient
    return final