from sage.all import *


def inverse(f, n):
    # calcul pour une serie formelle de la forme F=1+XG de son inverse modulo N
    if n == 1:
        return f.parent().one()
    else:
        G = inverse(f, (n + 1) // 2)
        return G + (1 - G.multiplication_trunc(f, n)).multiplication_trunc(G, n)


def log1(f, n):

    # calcul pour une serie formelle de la forme F=XG de log(1+F) modulo N
    D = f.derivative()
    I = inverse(1 + f, n)
    return (D.multiplication_trunc(I, n - 1)).integral()


def inverse2(f, n):

    # calcul l'inverse une serie formelle F avec F(0)!=0
    if n == 1:
        return f.parent()(1 / f(0))
    else:
        G = inverse(f, (n + 1) // 2)
        return G + (1 - G.multiplication_trunc(f, n)).multiplication_trunc(G, n)


def log2(f, n):
    a = f[0]
    return log(a) + log1((f / a - 1), n)


def exp(f, n):
    s = f.parent().one()
    for k in range(n):
        s += s.multiplication_trunc(f - log2(s, n), n)
    return s

def fsin(f,n):
    t=exp(f*CC.0,n)
    tt=exp(-f*CC.0,n)
    return -0.5*CC.0*(t-tt)

def fcos(f,n):
    t=exp(f*CC.0,n)
    tt=exp(-f*CC.0,n)
    return 0.5*(t+tt)

def ftan(f,n):
    return fsin(f,n)*inverse2(fcos(f,n),n)

def farctan(f,n):
    return 0.5*CC.0*log2(f.parent().one()-CC.0*f,n)-0.5*CC.0*log2(f.parent().one()+CC.0*f,n)
