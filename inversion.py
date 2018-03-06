from sage.all import *

R = PolynomialRing(QQ,'X') #Univariate Polynomial Ring in X over Rational Field
X = R.gen()


def inverse(F,N) :
    #calcul pour une serie formelle de la forme F=1+XG de son inverse modulo N
    if N == 1 :
        return 1
    else :
        G = inverse(F,(N +1)//2)
        return (G+(1-G.multiplication_trunc(F,N)).multiplication_trunc(G,N))

def Log1(F,N) :
    #calcul pour une serie formelle de la forme F=XG de log(1+F) modulo N
    D = F.derivative()
    I = inverse(1+F,N)
    return (D.multiplication_trunc(I,N-1)).integral()

