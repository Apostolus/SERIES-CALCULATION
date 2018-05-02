from sage.all import *


def inverse_un_series(f, n):
    # calcul pour une serie formelle de la forme F=1+XG de son inverse modulo N
    
    if n == 1:
        return f.parent().one()
    else:
        G = inverse_un_series(f, (n + 1) // 2)
        return G + (1 - G.multiplication_trunc(f, n)).multiplication_trunc(G, n)

def inverse_series(f, n):
    # calcul l'inverse une serie formelle F avec F(0)!=0
    
    if n == 1:
        return f.parent()(1 / f(0))
    else:
        G = inverse_series(f, (n + 1) // 2)
        return G + (1 - G.multiplication_trunc(f, n)).multiplication_trunc(G, n)

def log_un_series(f, n):
    # calcul pour une serie formelle de la forme F=XG de log(1+F) modulo N
    
    D = f.derivative()
    I = inverse_un_series(1 + f, n)
    return (D.multiplication_trunc(I, n - 1)).integral()

def log_zero_un_series(f, n):
    # calculer log(f) pour une serie formelle f, qui verifie forcement que f(0)=1

    return log_un_series((f - 1), n)

def log_zero_nonNull_series(f, n):
    # calculer log(f) pour une serie formelle f, qui verifie forcement que f(0)!=0
    a = f(0)
    return log(a) + log_un_series(f/a-1,n)

def exp_zero_zero_series(f, n):
    # calculer exp(f) pour une serie f qui suffit forcement que f(0)=0
    
    s = f.parent().one()
    for k in range(n):
        s += s.multiplication_trunc(f - log_zero_un_series(s, n), n)
    return s

def exp_series(f, n):
    a = f(0)
    return exp(a)*exp_zero_zero_series(f-a, n)

def fsin(f,n):
    t=exp_series(f*CC.0,n)
    tt=exp_series(-f*CC.0,n)
    return -0.5*CC.0*(t-tt)

def fcos(f,n):
    t=exp_series(f*CC.0,n)
    tt=exp_series(-f*CC.0,n)
    return 0.5*(t+tt)

def ftan(f,n):
    return fsin(f,n)*inverse_series(fcos(f,n),n)

def farctan(f,n):
    return 0.5*CC.0*log_zero_nonNull_series(f.parent().one()-CC.0*f,n)-0.5*CC.0*log_zero_nonNull_series(f.parent().one()+CC.0*f,n)
