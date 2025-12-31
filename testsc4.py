import scipy.integrate as integrate
import matplotlib.pyplot as plt
import itertools as iter
from sympy import bernoulli
import utils as ut
import conesfunc as con
import zetas as zet 
import dimtwo as dt
import resolplanecones as res
import cmath
import math

def testeq13(v0,v1,q):
    v = [v0,v1]
    Q= [1,q]
    final = zet.bernrec(v, Q)
    bq = bernoulli(q,0)
    if q%2 == 0:
        derived = bq*(v0**(1-q) - 1)
    if q%2 == 1:
        if q == 1:
            derived = 0
        else:
            derived = -bq
    print('realval = ', final)
    print('derived value =', derived)
    if derived == 0:
        print('DIVISON BY 0',  final,'----',derived)
    else:
        print('ratio =', final/derived) 
    return 0

def testteo41(v,q):
    for k in range(len(v)):
        finalZ = zet.Zperp(k, v, q)
        V = [v[j]  for j in range(len(v)) if j != k]
        Q = [q[j]  for j in range(len(q)) if j != k]
        derivedZ = zet.Zvq(V, Q)
        print('Z, k=', k , finalZ, derivedZ)
    finalY = 0
    derivedY = 0
    for k in range(0,len(v)):
        if q[k]==1:
            finalY += v[k]*zet.Yperp(k,v,q)
        if q[k]>1:
            Q = [q[j] if j != k else q[j]-1 for j in range(len(q))]
            derivedY += -v[k]*zet.Zvq(v, Q)
    print('Y',finalY, derivedY)
    return 0
 
def testremark421(v,q):
    final = 0
    Points = zet.vperp(v, 0)
    print(Points)
    for P in Points:
        quf = 0 
        for k in range(len(v)):
            qu = 1
            for j in range(len(v)):
                if j != k and P[j] != 0:
                    qu *= P[j]**(-q[j])
            quf += v[k]*qu  
        final += quf
    print('result = ', final)
    return 0

def testremark422(k,v,q):
    fs = [lambda x, s=qq: bernoulli(s, x) for qq in q]
    normqk = sum([q[j] for j in range(len(v)) if j != k])
    P = (2*math.pi*1j)**(normqk)
    for j in range(len(v)):
        if j != k :
            P = P/(math.factorial(q[j]))
    def prodf(x):
        prod = 1
        for j in range(len(v)):
            if j != k:
                prod *= fs[j](ut.frac(v[j]*x))
        return prod
    derived = ((-1)**(normqk))*P
    integral = integrate.quad(lambda x : prodf(x) , 0 ,1)
    derived *= integral[0]
    print(integral[0])
    final = zet.Zperp(k, v, q)
    print('realval = ', final)
    print('derived value =', derived)
    if derived == 0:
        print('DIVISON BY 0',  final,'----',derived)
    else:
        print('ratio =', final/derived)     
    return 0


def testprop43Y(k,q, N=50):
    r = len(q)-1
    x = 1
    for n in range(N):
        v = ut.randvs(r)
        cotaY, cotaZ, vl = ut.cotas(k, v, q)
        U = iter.product([-1,1],repeat = r+1)
        for u in U:
            #print(n)
            x = x + 1
            #print('x',x)
            y = zet.Yperpu(v,k,q,u)/cotaY
            plt.scatter(x,y)
        
    plt.show()
    return 1

def testprop43Z(k,q, N=50):
    r = len(q)-1
    x = 1
    for n in range(N):
        v = ut.randvs(r)
        cotaY, cotaZ , vl= ut.cotas(k, v, q)
        U = iter.product([-1,1],repeat = r+1)
        for u in U:
            #print(n)
            x = x + 1
            #print('x',x)
            y = zet.Zperpu(v,k,q,u)/cotaZ
            plt.scatter(x,y)
        
    plt.show()
    return 1

def testcor44(k,q,N=50):
    r = len(q)-1
    x = 1
    for n in range(N):
        v = ut.randvs(r)
        cotaY, cotaZ , vl= ut.cotas(k, v, q)
        cotaZeta = 2**(r+1) * (1+4*vl) * cotaZ 
        x = x + 1
        #print('x',x)
        y = zet.zetakvq(k, v, q)/cotaZeta
        plt.scatter(x,y)     
    plt.show()
    return 1



#v0 = 2
#v1 = 5
#q = 4
#testeq13(v0,v1,q)

v=[3,5,7]
q = [1,1,1]
#k = 1

#testteo41(v,q)

q = [1,1,1]
testremark421(v,q)

#testremark422(k,v,q)



