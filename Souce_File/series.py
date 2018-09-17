from sage.all import *

"""
#/usr/share/SageMath/src/sage/structure/element.pxd
from sage.structure.element import parent
#Nous ne connaissons pas encore le fichier qui contient la fonction "one" comme ".parent.one()"
#le fichier qui contient la fonction "zero" comme "parent.zero"
#/usr/share/SageMath/local/lib/python2.7/site-packages/IPython/core/autocall.py
from autocall import exit
#le fichier qui contient la fonction "multiplication_trunc"
#/usr/share/SageMath/local/lib/python2.7/site-packages/sage/calculus/functional.py
from calculus.functional import derivative
#/usr/share/SageMath/local/lib/python2.7/site-packages/sage/misc/functional.py
from misc.functional import integral
"""

"""
This package is built for the operations on series whose terme could be rational, real, complex and even terms of another ring. 
For example, x + 1/2*x^2 + 1/3*x^3, 1.0*x + 0.5*x^2 + 0.333*x^3, I*x + (0.5+0.5I)*x^2, and ax + 1/2*a^2*x^2, etc. 
And don't forget well declare the ring of each function before use them :
    for rational:
sage: Q = PolynomialRing(QQ, 'x')
sage: x = Q.gen()
    for real:
sage: R = PolynomialRing(RR, 'y')
sage: y = R.gen()
    for complex:
sage: C = PolynomialRing(CC, 'z')
sage: z = C.gen()
    for ring:
sage: AQ = PolynomialRing(Q, 'a')
sage: a = AQ.gen()
However, our goal and the principal advantages of this package is rapidly calculating the operations on rational, which means our functions can insure the result after operations on rational series is still rational, while accelerating calculions.
This package contains the basic operations on a serie: inverse, log, exp, arctan, tan, sin and cos.
In maths, generally, the number of terms of the exact result after these operations should be infinitied. So every time when we do these operations, we need give the level of precision, which is the max number of terms and presented by "n".
"""


def inverse_un_series(f, n):
    """
    This function returns the first "n" terms of the inverse of the serie "f" that satisfies f(0) is 1. In other word, it will return a serie g that satifies f * g = 1.
    If f(0) isn't 1, it will return raise an exception of error and exit.
    ALGORITHM:
        Newton iteration starting from (x^(-1))^(-1) = x
    
    EXAMPLE::
     
        sage: Q = PolynomialRing(QQ, 'x')
        sage: x = Q.gen()
    
        sage: P = x + 1
        sage: inverse_un_series(P, 10)
        -x^9 + x^8 - x^7 + x^6 - x^5 + x^4 - x^3 + x^2 - x + 1
    
    ::
    
        sage: P = x*x + 1
        sage: inverse_un_series(P, 10)
        x^8 - x^6 + x^4 - x^2 + 1
    ::
        sage: AQ = PolynomialRing(Q, 'a')
        sage: a = AQ.gen()
        
        sage: P = 1 + x*a + (x + 1)*a*a
        sage: inverse_un_series(P, 4)
        (-x^3 + 2*x^2 + 2*x)*a^3 + (x^2 - x - 1)*a^2 - x*a + 1
        
        sage: P.multiplication_trunc(inverse_un_series(P, 4), 4)
        1
    """
    
    if(not f(0).is_one()):
        raise ZeroDvisionError("Error: f(0) must be equal to 1")

    if n == 1:
        return f.parent().one()
    else:
        G = inverse_un_series(f, (n + 1) // 2)
        return G + (1 - G.multiplication_trunc(f, n)).multiplication_trunc(G, n)

def inverse_series(f, n):
	"""
	Calculates the inverse of a serie F with F(0) invertible	
	"""
    if(f(0).is_one()):
        return inverse_un_series(f,n)

    if(f(0)==0):
        raise ZeroDivisionError("Error: f(0) is 0")

    if (not f(0).is_unit()):
        raise ZeroDivisionError("Error: f(0) is not invertible")

    if n==1:
        return f.parent()(f(0).inverse_of_unit())
    else:
        G = inverse_series(f, (n + 1)//2)
        return G + (1 - G.multiplication_trunc(f, n)).multiplication_trunc(G, n)

def log_un_series(f, n):
	"""
	This function is the heart part of the next function to calculate logarithm of a polynom, it executes the algorithm which uses the formula log(x) = integral(1/x). But, it only accepts the polynome in form of F=XG, G is a polynome. Finally, it returns log(f) in form of a serie.
	"""
    if(f(0)!=0):
        print("Error: f(0) must be equal to 0!")
        sys.exit(1)

    D = f.derivative()
    I = inverse_un_series(1 + f, n)
    return (D.multiplication_trunc(I, n - 1)).integral()

def log_zero_un_series(f, n):
    """
    This function returns the firsr "n" termes of the logarithm of serie "f" that satifies f(0) = 1. In other word, it will return a serie g that satifies e^g = f.
    If f(0) isn't 1, it will raise an exception and exit.
    ALGORITHM:
        Uses the formula log(x) = integral(1/x)
    EXAMPLE::
    
        sage: Q = PolynomialRing(QQ, 'x')
        sage: x = Q.gen()
        sage: P = x + 1
        sage: log_zero_un_series(P, 10)
        1/9*x^9 - 1/8*x^8 + 1/7*x^7 - 1/6*x^6 + 1/5*x^5 - 1/4*x^4 + 1/3*x^3 - 1/2*x^2 + x
    ::
        sage: P = x*x + 1
        sage: log_zero_un_series(P, 10)
        -1/4*x^8 + 1/3*x^6 - 1/2*x^4 + x^2
    
    ::
        sage: AQ = PolynomialRing(Q, 'a')
        sage: a = AQ.gen()
        sage: P = a*x + 1
        sage: log_zero_un_series(P, 10)
        1/9*x^9*a^9 - 1/8*x^8*a^8 + 1/7*x^7*a^7 - 1/6*x^6*a^6 + 1/5*x^5*a^5 - 1/4*x^4*a^4 + 1/3*x^3*a^3 - 1/2*x^2*a^2 + x*a
    """
    
    if(f(0)!=f.parent().one()):
        raise ZeroDivisionError("Error: f(0) must be equal to 1")

    return log_un_series((f - 1), n)

def log_zero_nonNull_series(f, n):
	"""
	This function is a more general logarithm function, which means there is only one mathematical limit on the polynome f, which is f(0) != 0, instead of the need of "f(0) = 1" in the last function "log_zero_un_series(f, n)". It can be used for any polynome f, even if f(0) != 1, for example a, because we can separe this operation in two parts: log(a) and log(f - a), finally it will return log(a) + log(f - a). But, please pay your attention, this function cannot insure the result is still on the same ring of f, because of log(a). So, the advise is to use it when f bases on real or complex.
	"""
    if(f(0)==0):
        raise ZeroDivisionError("Error: log(0) is not defined")
    
    a = f(0)
    return log(a) + log_un_series(f/a - 1, n)

def exp_zero_zero_series(f, n):
    """
    This function returns the first "n" termes of the exponent of serie "f" that satisfies f(0) = 0. In other words, it will return a serie g that satisfies e^f = g.
    If f(0) isn't 0, it will raise an exception and exit.
    ALGORITHM:
        Newton iteration starting from log(exp(x)) = x.
    EXAMPLE::
        sage: Q = PolynomialRing(QQ, 'x')
        sage: x = Q.gen()
        sage: P = x
        sage: exp_zero_zero_series(P, 10)
        1/362880*x^9 + 1/40320*x^8 + 1/5040*x^7 + 1/720*x^6 + 1/120*x^5 + 1/24*x^4 + 1/6*x^3 + 1/2*x^2 + x + 1
    ::
        P = x + 2*x*x + 3*x*x*x
        sage: exp_zero_zero_series(P, 10)
        1626437/72576*x^9 + 794081/40320*x^8 + 15431/1008*x^7 + 9661/720*x^6 + 1181/120*x^5 + 145/24*x^4 + 31/6*x^3 + 5/2*x^2 + x + 1
    
    ::
        sage: AQ = PolynomialRing(Q, 'a')
        sage: a = AQ.gen()
        sage: P = a*x
        sage: exp_zero_zero_series(P, 10)
        1/362880*x^9*a^9 + 1/40320*x^8*a^8 + 1/5040*x^7*a^7 + 1/720*x^6*a^6 + 1/120*x^5*a^5 + 1/24*x^4*a^4 + 1/6*x^3*a^3 + 1/2*x^2*a^2 + x*a + 1
    
    ::
        P = x^2*a^2 + x*a
        sage: log_zero_un_series(exp_zero_zero_series(P, 10), 10)==P
        True
    """
    if(f(0)!=f.parent().zero()):
        print("Error: f(0) must be equal to 0!")
        sys.exit(1)
    
    s = f.parent().one()
    #optimisation sur le nombre de fois de l'iterations
    #r = int(log(n, 2))
    #r=n
    
    #for k in range(r):
    #    s += s.multiplication_trunc(f - log_zero_un_series(s, n), n)
    deg = 1
    while(deg<n) :
        s += s.multiplication_trunc(f - log_zero_un_series(s, deg + 1), deg)
        deg = deg*2
    s += s.multiplication_trunc(f - log_zero_un_series(s, deg+1), n)
    return s

def exp_series(f, n):
	"""
	This function is a more general exponent function, which means there is no limit on the polynome f, instead of the need of "f(0) = 0" in the last function "exp_zero_zero_series(f, n)". It can be used for any polynome f, even if f(0) != 0, for example a, because we can separe this operation in two parts: exp(a) and exp(f - a), finally it will return exp(a)*exp(f - a). But, please pay your attention, this function cannot insure the result is still on the same ring of f, because of exp(a). So, the advise is to use it when f bases on real or complex.
	"""

    a = f(0)
    return exp(a)*exp_zero_zero_series(f - a, n)

def arctan_series(f, n):
    """
    This function returns the first "n" termes of arctan of serie f that satisfies f(0) = 0.
    If f(0) isn't 0, it will raise an exception and exit.
    ALGORITHM:
        Uses the formula arctan(x) = integral(1/(1 + x))
    EXAMPLE::
    
        sage: Q = PolynomialRing(QQ, 'x')
        sage: x = Q.gen()
        sage: P = x
        sage: arctan_series(P, 10)
        1/9*x^9 - 1/7*x^7 + 1/5*x^5 - 1/3*x^3 + x
    ::
        sage: P = x + 2*x*x
        sage: arctan_series(P,10) == P._atan_series(10)
        True
    ::
        sage: AQ = PolynomialRing(Q, 'a')
        sage: a = AQ.gen()
        
        sage: P = x*x*a*a + 2*x*a
        sage: arctan_series(P, 10)
        -334/9*x^9*a^9 - 56*x^8*a^8 - 16/7*x^7*a^7 + 47/3*x^6*a^6 + 22/5*x^5*a^5 - 4*x^4*a^4 - 8/3*x^3*a^3 + x^2*a^2 + 2*x*a
    """

    g = f.derivative()
    h = f.parent().one() + f.multiplication_trunc(f, n)
    t = inverse_series(h, n)
    r = g.multiplication_trunc(t, n-1)
    return r.integral()

def tan_series(f, n):
    """
    This function returns the first "n" termes of tan of serie f that satisfies f(0) = 0.
    If f(0) isn't 0, it will raise an exception and exit.
    ALGORITHM:
        Newton iteration starting from arctan(tan(x)) = x
    EXAMPLE::
        sage: Q=PolynomialRing(QQ,'x')
        sage: x=Q.gen()
    
        sage: P=x
        sage: tan_series(P,10)
        62/2835*x^9 + 17/315*x^7 + 2/15*x^5 + 1/3*x^3 + x
    ::
        sage: P=x^2 + 2*x
        sage: tan_series(P,10)
        27668/567*x^9 + 1328/45*x^8 + 5536/315*x^7 + 11*x^6 + 94/15*x^5 + 4*x^4 + 8/3*x^3 + x^2 + 2*x
    ::
        sage: AQ=PolynomialRing(Q,'a')
        sage: a=AQ.gen()
        sage: P=x*x*a*a+2*x*a
        sage: tan_series(P,10)
        27668/567*x^9*a^9 + 1328/45*x^8*a^8 + 5536/315*x^7*a^7 + 11*x^6*a^6 + 94/15*x^5*a^5 + 4*x^4*a^4 + 8/3*x^3*a^3 + x^2*a^2 + 2*x*a
    """
    if(f(0)!=0):
        print("Error: f(0) must be equal to 0!")
        sys.exit(1)

    un = f.parent().one()
    s = f
    deg = 1
    while(deg<n):
        temps1 = f - arctan_series(s, deg+1)
        temps2 = un + s.multiplication_trunc(s, deg+1)
        s += temps1.multiplication_trunc(temps2, deg)
        deg*=2
    temps1 = f - arctan_series(s, deg+1)
    temps2 = un + s.multiplication_trunc(s, n)
    s += temps1.multiplication_trunc(temps2, n)
    return s

def sin_series(f, n):
    """
    This function returns the first "n" termes of the sinus of serie f that satisfies f(0) = 0.
    If f(0) isn't 0, it will raise an exception and exit.
    ALGORITHM:
        Uses the formula sin(x) = 2*tan(x/2)/(1 + tan(x/2))
    EXAMPLE::
        sage: Q = PolynomialRing(QQ,'x')
        sage: x = Q.gen()
    
        sage: P = x
        sage: sin_series(P, 10)
        1/362880*x^9 - 1/5040*x^7 + 1/120*x^5 - 1/6*x^3 + x
    ::
        sage: P = x*x + 2*x
        sage: sin_series(P, 10)
        -551/11340*x^9 + 11/45*x^8 + 202/315*x^7 + 1/2*x^6 - 11/15*x^5 - 2*x^4 - 4/3*x^3 + x^2 + 2*x
    ::
        sage: AQ = PolynomialRing(Q, 'a')
        sage: a = AQ.gen()
        
        sage: H = x*x*a*a + 2*a*x
        sage: sin_series(H, 10)
        -551/11340*x^9*a^9 + 11/45*x^8*a^8 + 202/315*x^7*a^7 + 1/2*x^6*a^6 - 11/15*x^5*a^5 - 2*x^4*a^4 - 4/3*x^3*a^3 + x^2*a^2 + 2*x*a
    """

    temps = tan_series(f/2, n)
    temps2 = f.parent().one() + temps.multiplication_trunc(temps, n)
    return 2*temps.multiplication_trunc(inverse_series(temps2, n), n)

def cos_series(f, n):
    """
    This function returns the first "n" termes of the cosinus of serie f that satisfies f(0) = 0.
    If f(0) isn't 0, it will raise an exception and exit.
    
    ALGORITHM:
        Uses the formula cos(x) = (1 - tan(x/2)^2)/(1 + tan(x/2)^2)
		
    EXAMPLE::
        
        sage: Q = PolynomialRing(QQ, 'x')
        sage: x = Q.gen()
        
        sage: P = x
        sage: cos_series(P, 10)
        1/40320*x^8 - 1/720*x^6 + 1/24*x^4 - 1/2*x^2 + 1
    ::
        sage: P = x*x + 2*x
        sage: cos_series(P, 10)
        -62/315*x^9 - 719/2520*x^8 + 1/15*x^7 + 41/45*x^6 + 4/3*x^5 + 1/6*x^4 - 2*x^3 - 2*x^2 + 1
    ::
        sage: AQ = PolynomialRing(Q, 'a')
        sage: a = AQ.gen()
        sage: P = a*a*x*x + 2*a*x
        sage: cos_series(P, 10)
        -62/315*x^9*a^9 - 719/2520*x^8*a^8 + 1/15*x^7*a^7 + 41/45*x^6*a^6 + 4/3*x^5*a^5 + 1/6*x^4*a^4 - 2*x^3*a^3 - 2*x^2*a^2 + 1
    """

    un = f.parent().one()
    temps = tan_series(f/2, n)
    temps = temps.multiplication_trunc(temps, n)
    temps2 = un - temps
    temps3 = un + temps
    return temps2.multiplication_trunc(inverse_series(temps3, n), n)
