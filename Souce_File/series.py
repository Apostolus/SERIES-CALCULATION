from sage.all import *

"""
sage: Q=PolynomialRing(QQ,'x')
sage: x=Q.gen()
sage: R=PolynomialRing(RR,'y')
sage: y=R.gen()
sage: C=PolynomialRing(CC,'z')
sage: z=C.gen()
sage: AQ=PolynomialRing(Q,'a')
sage: a=AQ.gen()
"""


def inverse_un_series(f, n):
    # calcul pour une serie formelle de la forme F=1+XG de son inverse modulo N
    """
    EXAMPLE::
        sage: P=x+1
        sage: inverse_un_series(P,10)
        -x^9 + x^8 - x^7 + x^6 - x^5 + x^4 - x^3 + x^2 - x + 1

        sage: P=x*x+1
        sage: inverse_un_series(P,10)
        x^8 - x^6 + x^4 - x^2 + 1

        sage: AQ
        Univariate Polynomial Ring in a over Univariate Polynomial Ring in x over Rational Field

        sage: P=1+x*a+(x+1)*a*a
        sage: inverse_un_series(P,4)
        (-x^3 + 2*x^2 + 2*x)*a^3 + (x^2 - x - 1)*a^2 - x*a + 1
        sage: P.multiplication_trunc( inverse_un_series(P,4),4)
        1

    """
    
    if(f(0)!=f.parent().one()):
        print("Erreur! f(x) n'est pas sous forme de f=1+xG!")
        sys.exit(1)

    if n == 1:
        return f.parent().one()
    else:
        G = inverse_un_series(f, (n + 1) // 2)
        return G + (1 - G.multiplication_trunc(f, n)).multiplication_trunc(G, n)

def inverse_series(f, n):
    # calcul l'inverse une serie formelle F avec F(0)!=0
    if(f(0)==f.parent().one()):
        return inverse_un_series(f,n)
    if(f(0)==0):
        print("Erreur! f(0) ne doit pas etre zero!")
        sys.exit(1)

    if n == 1:
        return f.parent()(1 / f(0))
    else:
        G = inverse_series(f, (n + 1) // 2)
        return G + (1 - G.multiplication_trunc(f, n)).multiplication_trunc(G, n)

def log_un_series(f, n):
    # calcul pour une serie formelle de la forme F=XG de log(1+F) modulo N
    if(f(0)!=0):
        print("Erreur! f(0) doit etre 0!")
        sys.exit(1)

    D = f.derivative()
    I = inverse_un_series(1 + f, n)
    return (D.multiplication_trunc(I, n - 1)).integral()

def log_zero_un_series(f, n):
    # calculer log(f) pour une serie formelle f, qui verifie forcement que f(0)=1
    """
    EXAMPLE::
        sage: P=x+1
        sage: log_zero_un_series(P,10)
        1/9*x^9 - 1/8*x^8 + 1/7*x^7 - 1/6*x^6 + 1/5*x^5 - 1/4*x^4 + 1/3*x^3 - 1/2*x^2 + x

        sage: P=x*x+1
        sage: log_zero_un_series(P,10)
        -1/4*x^8 + 1/3*x^6 - 1/2*x^4 + x^2
        sage: log_zero_un_series(P,10)==P._log_series(10)
        True

        sage: P=a*x+1
        sage: log_zero_un_series(P,10)
        1/9*x^9*a^9 - 1/8*x^8*a^8 + 1/7*x^7*a^7 - 1/6*x^6*a^6 + 1/5*x^5*a^5 - 1/4*x^4*a^4 + 1/3*x^3*a^3 - 1/2*x^2*a^2 + x*a
    """
    
    if(f(0)!=f.parent().one()):
        print("Erreur! f(0) doit etre 1!")
        sys.exit(1)

    return log_un_series((f - 1), n)

def log_zero_nonNull_series(f, n):
    # calculer log(f) pour une serie formelle f, qui verifie forcement que f(0)!=0
    if(f(0)==0):
        print("Erreur! log(0) n'est pas defini!")
        sys.exit(1)
    
    a = f(0)
    return log(a) + log_un_series(f/a-1,n)

def exp_zero_zero_series(f, n):
    # calculer exp(f) pour une serie f qui suffit forcement que f(0)=0
    """
    EXAMPLE::
        sage: P=x
        sage: exp_zero_zero_series(P,10)
        1/362880*x^9 + 1/40320*x^8 + 1/5040*x^7 + 1/720*x^6 + 1/120*x^5 + 1/24*x^4 + 1/6*x^3 + 1/2*x^2 + x + 1

        P=x+2*x*x+3*x*x*x
        sage: exp_zero_zero_series(P,10)==P._exp_series(10)
        True

        sage: P=a*x
        sage: exp_zero_zero_series(P,10)
        1/362880*x^9*a^9 + 1/40320*x^8*a^8 + 1/5040*x^7*a^7 + 1/720*x^6*a^6 + 1/120*x^5*a^5 + 1/24*x^4*a^4 + 1/6*x^3*a^3 + 1/2*x^2*a^2 + x*a + 1
        
        P=x^2*a^2 + x*a
        sage: log_zero_un_series(exp_zero_zero_series(P,10),10)==P
        True
    """
    if(f(0)!=0):
        print("Erreur! f(0) doit etre 0!")
        sys.exit(1)
    
    s = f.parent().one()
    #optimisation sur le nombre de fois de l'iterations
    #r = int(log(n,2))
    r=n
    for k in range(r):
        s += s.multiplication_trunc(f - log_zero_un_series(s, n), n)
    return s

def exp_series(f, n):
    # calculer exp(f) pour une serie f

    a = f(0)
    return exp(a)*exp_zero_zero_series(f-a, n)

def arctan_series(f,n):
    # calculer arctan(f) pour une serie f
    """
    EXAMPLE::
        sage: P=x
        sage: arctan_series(P,10)
        1/9*x^9 - 1/7*x^7 + 1/5*x^5 - 1/3*x^3 + x

        sage: P=x+2*x*x
        sage: arctan_series(P,10)==P._atan_series(10)
        True

        sage: arctan_series(P,10)
        -334/9*x^9*a^9 - 56*x^8*a^8 - 16/7*x^7*a^7 + 47/3*x^6*a^6 + 22/5*x^5*a^5 - 4*x^4*a^4 - 8/3*x^3*a^3 + x^2*a^2 + 2*x*a
    """

    g=f.derivative()
    h=f.parent().one()+f.multiplication_trunc(f,n)
    t=inverse_series(h,n)
    r=g.multiplication_trunc(t,n-1)
    return r.integral()

def tan_series(f,n):
    #calculer tan(f) pour une serie f
    """
    EXAMPLE::
        sage: P=x
        sage: tan_series(P,10)
        62/2835*x^9 + 17/315*x^7 + 2/15*x^5 + 1/3*x^3 + x

        sage: P=x^2 + 2*x
        sage: tan_series(P,10)==P._tan_series(10)
        True
        
        sage: tan_series(P,10)
        27668/567*x^9*a^9 + 1328/45*x^8*a^8 + 5536/315*x^7*a^7 + 11*x^6*a^6 + 94/15*x^5*a^5 + 4*x^4*a^4 + 8/3*x^3*a^3 + x^2*a^2 + 2*x*a
    """
    if(f(0)!=0):
        print("Erreur! f(0) doit etre 0!")
        sys.exit(1)

    un = f.parent().one()
    s = f
    for k in range(n):
        temps1 = f - arctan_series(s,n)
        temps2 = un + s.multiplication_trunc(s,n)
        s += temps1.multiplication_trunc(temps2,n)
    return s

def sin_series(f,n):
    # calculer sin(f) pour une serie f
    """
    EXAMPLE::
        sage: P=x
        sage: sin_series(P,10)
        1/362880*x^9 - 1/5040*x^7 + 1/120*x^5 - 1/6*x^3 + x

        sage: P=x*x+2*x
        sage: sin_series(P,10)==P._sin_series(10)
        True
        
        sage: H=x*x*a*a+2*a*x
        sage: sin_series(H,10)
        -551/11340*x^9*a^9 + 11/45*x^8*a^8 + 202/315*x^7*a^7 + 1/2*x^6*a^6 - 11/15*x^5*a^5 - 2*x^4*a^4 - 4/3*x^3*a^3 + x^2*a^2 + 2*x*a

    """

    temps = tan_series(f/2,n)
    temps2 = f.parent().one() + temps.multiplication_trunc(temps,n)
    return 2*temps.multiplication_trunc(inverse_series(temps2,n),n)

def cos_series(f,n):
    # calculer cos(f) pour une serie f
    """
    EXAMPLE::
        sage: P=x
        sage: cos_series(P,10)
        1/40320*x^8 - 1/720*x^6 + 1/24*x^4 - 1/2*x^2 + 1

        sage: P=x*x+2*x
        sage: cos_series(P,10)
        -62/315*x^9 - 719/2520*x^8 + 1/15*x^7 + 41/45*x^6 + 4/3*x^5 + 1/6*x^4 - 2*x^3 - 2*x^2 + 1
        sage: cos_series(P,10)==P._cos_series(10)
        True

        sage: P=a*a*x*x+2*a*x
        sage: cos_series(P,10)
        -62/315*x^9*a^9 - 719/2520*x^8*a^8 + 1/15*x^7*a^7 + 41/45*x^6*a^6 + 4/3*x^5*a^5 + 1/6*x^4*a^4 - 2*x^3*a^3 - 2*x^2*a^2 + 1

    """

    un = f.parent().one()
    temps = tan_series(f/2,n)
    temps = temps.multiplication_trunc(temps,n)
    temps2 = un - temps
    temps3 = un + temps
    return temps2.multiplication_trunc(inverse_series(temps3,n),n)

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
