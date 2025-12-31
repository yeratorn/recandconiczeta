import resolplanecones as res
import conesfunc as con
from sympy import bernoulli
import math


def sigma(k,l,q):
    if k>len(q) or l>len(q):
        return 'dimension error'
    U = [[1,-1,-1],[1,-1,1],[1,1,-1]]
    final = 1
    for j in range(len(q)):
        if j != k:
            final *= U[l][j]**(q[j])
    return final

def Zvalued(v,q,k):
    v0,v1,v2 = v
    Vect = [[0,v2,v1],[v2,0,v0],[v1,v0,0]]
    ran = [0,1,2]
    sig = 0
    for l in ran:
        if l != k:
            sig += sigma(k,l,q)
    quoti = 1 
    for j in [0,1,2]:
        if j != k:
            quoti *= Vect[k][j]**(-q[j])
    derivedval = quoti*sig
    return derivedval
  

def cones(v):
    return [[[v[2],0,v[0]],[v[1],v[0],0]],[[v[1],v[0],0],[0,v[2],v[1]]],[[0,v[2],v[1]],[v[2],0,v[0]]]]

def Qvaluekl(v,q,k,l):
    CONES = cones(v)
    V1 = CONES[l][0]
    V2 = CONES[l][1]
    RES = res.res2dcone(V1,V2)
    s = len(RES)-1
    suma1 = 0
    Qk = 0 
    for i in range(1,s):
        quoti = 1
        for j in [0,1,2]:
            if j != k:
                quoti *= RES[i][j]**(-q[j])
        suma1 += quoti
    suma2 = 0
    qk = [q[j] if j != k else 0 for j in [0,1,2]]
    for i in range(0,s):
        suma2 += con.zetaconegcd([RES[i],RES[i+1]],qk)
    Qk += suma1+suma2
    return Qk

def delta2d(q):
    unos = 0
    prodb = 1
    for j in range(len(q)):
        if q[j] == 1:
            unos += 1 
        if q[j] != 1:
            prodb *= bernoulli(q[j],0)
    un = 1 
    un = (1 - (-1)**unos)*(2**(-unos))
    return un*prodb

def constant(q,k):
    q[k] = 0
    C = -bernoulli(sum(q),0)
    C *= math.factorial(sum(q))**(-1)
    for qj in q:
        C *= math.factorial(qj)
    return C

def constquno(q):
    suma = sum([q[j] for j in range(len(q))])
    suma = suma - 1
    prod = -bernoulli(suma, 0)/(math.factorial(suma))
    for qj in q:
        if qj>1:
            prod *= math.factorial(qj)
    return prod

def Qvalued(v,q,k):
    Qval = 0
    for l in [0,1,2]:
        Qval += sigma(k,l,q)*Qvaluekl(v,q,k,l)
    return Qval