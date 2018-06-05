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
    #optimisation sur le nombre de fois de l'iterations
    r = int(log(n,2))
    for k in range(r):
        s += s.multiplication_trunc(f - log_zero_un_series(s, n), n)
    return s

def exp_series(f, n):
    a = f(0)
    return exp(a)*exp_zero_zero_series(f-a, n)

def arctan_series(f,n):
    # calculer arctan(f) pour une serie f

    g=f.derivative()
    h=f.parent().one()+f.multiplication_trunc(f,n)
    t=inverse_series(h,n)
    r=g.multiplication_trunc(t,n-1)
    return r.integral()

def tan_series(f,n):
    #calculer tan(f) pour une serie f
    
    un = f.parent().one()
    s = f
    for k in range(n):
        temps1 = f - arctan_series(s,n)
        temps2 = un + s.multiplication_trunc(s,n)
        s += temps1.multiplication_trunc(temps2,n)
    return s

def sin_series(f,n):
    # calculer sin(f) pour une serie f

    temps = tan_series(f/2,n)
    temps2 = f.parent().one() + temps.multiplication_trunc(temps,n)
    return 2*temps.multiplication_trunc(inverse_series(temps2,n),n)

def cos_series(f,n):
    # calculer cos(f) pour une serie f

    un = f.parent().one()
    temps = tan_series(f/2,n)
    temps = temps.multiplication_trunc(temps,n)
    return (un + temps).multiplication_trunc(inverse_series(un - temps,n),n)

#--------------------------------------------------#
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
