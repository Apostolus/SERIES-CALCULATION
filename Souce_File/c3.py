from sage.all import *
from series import *
from random import *
from time import *
from math import *
import matplotlib.pyplot as plt


Q=PolynomialRing(QQ,'x')
x=Q.gen()
R=PolynomialRing(RR,'y')
y=R.gen()
C=PolynomialRing(CC,'z')
z=C.gen()
AQ=PolynomialRing(Q,'a')
a=AQ.gen()


def poly_random(switch,nbTermes,maxInt):
    """
    C'est une fonction qui revoie un polynome random avec des restrictions donnees:
    switch = le code de mode: 1->rationnel, 2->reel, 3->imaginaire, 4->anneau de l'anneau quotient
    nbTermes = le nombre de termes dans un polynome
    maxInt = la valeur absolue maximale de numerateur et denominateur pour le mode rationnel
             ou celle de coefficiences pour le mode reel
             ou celle de parties reelle et imaginaire de coefficience pour le mode imaginaire
    """
    P = 0

    if(switch == 1):
        for k in range(nbTermes):
            P = P*x + (int(2*maxInt*random())-maxInt)/(int(maxInt*random())+1)
    elif(switch == 2):
        for k in range(nbTermes):
            P = P*y + 2*maxInt*random()-maxInt
    elif switch == 3: 
        for k in range(nbTermes):
            P = P*z + 2*maxInt*random()-maxInt + (2*maxInt*random()-maxInt)*i
    else:
        print("Mode non-connu")
        sys.exit(1)

    return P


def polyCorrege(polynome, type):
    """
    type = 0: f(0)=0
    type = 1: f(0)=1
    """
    res = polynome
    if type==0:
        res = res - res(0) + polynome.parent().zero()
        return res
    if type==1:
        res = res - res(0) + polynome.parent().one()
        return res
    return polynome

#moyen = tempsMoyen(nbRepete,fun,switch,k,maxInt,varFixee)	
#poly_random(switch,nbTermes,maxInt)

def tempsMoyen(nb,fun,type,switch,nbTermes,maxInt,precision):
    """
    C'est une fonction qui teste le temps moyen d'effectuer la fonction 'fun(*paras)' avec 'nb' fois des iterations
    """
    moyen = 0
    for k in range(nb):
        p = poly_random(switch,nbTermes,maxInt)
        p = polyCorrege(p,type)
        debut = time()
        fun(p,precision)
        fin = time()
        moyen = moyen + fin - debut

    return moyen/nb


def reportWriter(path,mode,varFixee,nbIteration,nbRepete,fun,type,switch,maxInt):
    """
    C'est une fonction qui ecrit le resultat du test d'une fonction dans un fichier "path". Elle utilise la fonction tempsMoyen, et donc elle a besoin de paramettres de cette fonction, c'est-a-dire (nb,fun,*paras) de tempsMoyen
    Pour tester une fonction on a besoin de fixer la variable, donc on pose que mode 1 est de tester l'influence de nombre de termes donc on fixe le niveau de precision "varFixee"; et on pose que mode 2est de tester l'influen de niveau de precision, donc on fixe le nombre de termes.
    """
    f = open(path,"w")
    if(mode == 1):
        for k in range(1,nbIteration+1):
            print("\t\t%d termes, precision fixee = %d, on effectuera %d fois des iterations" %(k,varFixee,nbRepete))
            moyen = tempsMoyen(nbRepete,fun,type,switch,k,maxInt,varFixee)
            f.write("%d %.15f\n" %(k,moyen))
            print("\t\t  Temps moyen = %.15f seconds" %(moyen))

    elif(mode == 2):
        for k in range(2,nbIteration+2):
            print("\t\t%d fois de precision, nombre de termes fixe = %d, on effectuera %d fois des iterations" %(k,varFixee,nbRepete))
            moyen = tempsMoyen(nbRepete,fun,type,switch,varFixee,maxInt,k)
            f.write("%d %.15f\n" %(k,moyen))
            print("\t\t  Temps moyen = %.15f seconds" %(moyen))

    elif(mode == 3):
        for k in range(1,nbIteration):
            print("\t\t%d termes, precision fixee = %d, on effectuera %d fois des iterations" %(2**k,varFixee,nbRepete))
            moyen = tempsMoyen(nbRepete,fun,type,switch,2**k,maxInt,varFixee)
            f.write("%d %.15f\n" %(k,moyen))
            print("\t\t  Temps moyen = %.15f seconds" %(moyen))
    
    elif(mode == 4):
        for k in range(1,nbIteration):
            print("\t\t%d fois de precision, nombre de termes fixe = %d, on effectuera %d fois des iterations" %(2**k,varFixee,nbRepete))
            moyen = tempsMoyen(nbRepete,fun,type,switch,varFixee,maxInt,2**k)
            f.write("%d %.15f\n" %(k,moyen))
            print("\t\t  Temps moyen = %.15f seconds" %(moyen))
    else:
        print("Erreur: mode inconnu!")
        f.close()
        sys.exit(1)

    f.close()
    print("\n")
	

def test(nbIteration, switch, num, shift):
    """
    C'est une fonction qui teste le comportement de toutes les fonctions demandees dans le fichier series.py
    nbIteration c'est le nombre de fois que 
    """
    rootpath = "/usr/share/SageMath/2I013/data/"
    maxInt = 100
    #nbIteration = 100
    nbRepete = 3
    moyen = 0
    file_abs = ""
    nbTermesFixe = 5
    niveauPrecisionFixe = 5
    type = 0

    print("num=%d\n" %(num))
    if(num>=1000000):
        #---------------------------------inverse_un------------------------------------#
        type = 1
        print("Test de inverse_un_series, maxInt=%d\n" %(maxInt))
        ##----------------variable = nombre de termes-------------------##
        print("\tinverse_un TERMES")
        file_abs = rootpath + "inverse_un_terme.txt"
        reportWriter(file_abs,1+shift,niveauPrecisionFixe,nbIteration,nbRepete,inverse_un_series,type,switch,maxInt)
        ##---------------variable = niveau de precision-----------------##
        print("\tinverse_un PRECISION")
        file_abs = rootpath + "inverse_un_precision.txt"
        reportWriter(file_abs,2+shift,nbTermesFixe,nbIteration,nbRepete,inverse_un_series,type,switch,maxInt)
        
        num = num - 1000000

    print("num=%d\n" %(num))
    if(num>=100000):
        #------------------------------log_zero_un_series---------------------------------#
        type = 1
        print("Test de log_zero_un_series, maxInt=%d\n" %(maxInt))    
        ##----------------variable = nombre de termes-------------------##
        print("\tlog_zero_un TERMES")
        file_abs = rootpath + "log_zero_un_terme.txt"
        reportWriter(file_abs,1+shift,niveauPrecisionFixe,nbIteration,nbRepete,log_zero_un_series,type,switch,maxInt)
        ##---------------variable = niveau de precision-----------------##
        print("\tlog_zero_un PRECISION")
        file_abs = rootpath + "log_zero_un_precision.txt"
        reportWriter(file_abs,2+shift,nbTermesFixe,nbIteration,nbRepete,log_zero_un_series,type,switch,maxInt)
        
        num = num - 100000

    print("num=%d\n" %(num))
    if(num>=10000):
        #-----------------------------exp_zero_zero_series--------------------------------#
        type = 0
        print("Test de exp_zero_zero_series, maxInt=%d\n" %(maxInt))    
        ##----------------variable = nombre de termes-------------------##
        print("\texp_zero_zero TERMES")
        file_abs = rootpath + "exp_zero_zero_terme.txt"
        reportWriter(file_abs,1+shift,niveauPrecisionFixe,nbIteration,nbRepete,exp_zero_zero_series,type,switch,maxInt)
        ##---------------variable = niveau de precision-----------------##
        print("\texp_zero_zero PRECISION")
        file_abs = rootpath + "exp_zero_zero_precision.txt"
        reportWriter(file_abs,2+shift,nbTermesFixe,nbIteration,nbRepete,exp_zero_zero_series,type,switch,maxInt)
        
        num = num - 10000
    
    print("num=%d\n" %(num))
    if(num>=1000):
        #--------------------------------arctan_series-----------------------------------#
        type = 0
        print("Test de arctan_series, maxInt=%d\n" %(maxInt))    
        ##----------------variable = nombre de termes-------------------##
        print("\tarctan TERMES")
        file_abs = rootpath + "arctan_terme.txt"
        reportWriter(file_abs,1+shift,niveauPrecisionFixe,nbIteration,nbRepete,arctan_series,type,switch,maxInt)
        ##---------------variable = niveau de precision-----------------##
        print("\tarctan PRECISION")
        file_abs = rootpath + "arctan_precision.txt"
        reportWriter(file_abs,2+shift,nbTermesFixe,nbIteration,nbRepete,arctan_series,type,switch,maxInt)

        num = num - 1000

    print("num=%d\n" %(num))
    if(num>=100):
        #----------------------------------tan_series-------------------------------------#
        type = 0
        print("Test de tan_series, maxInt=%d\n" %(maxInt))    
        ##----------------variable = nombre de termes-------------------##
        print("\ttan TERMES")
        file_abs = rootpath + "tan_terme.txt"
        reportWriter(file_abs,1+shift,niveauPrecisionFixe,nbIteration,nbRepete,tan_series,type,switch,maxInt)
        ##---------------variable = niveau de precision-----------------##
        print("\ttan PRECISION")
        file_abs = rootpath + "tan_precision.txt"
        reportWriter(file_abs,2+shift,nbTermesFixe,nbIteration,nbRepete,tan_series,type,switch,maxInt)

        num = num -100

    print("num=%d\n" %(num))
    if(num>=10):
        #----------------------------------sin_series-------------------------------------#
        type = 0
        print("Test de sin_series, maxInt=%d\n" %(maxInt))    
        ##----------------variable = nombre de termes-------------------##
        print("\tsin TERMES")
        file_abs = rootpath + "sin_terme.txt"
        reportWriter(file_abs,1+shift,niveauPrecisionFixe,nbIteration,nbRepete,sin_series,type,switch,maxInt)
        ##---------------variable = niveau de precision-----------------##
        print("\tsin PRECISION")
        file_abs = rootpath + "sin_precision.txt"
        reportWriter(file_abs,2+shift,nbTermesFixe,nbIteration,nbRepete,sin_series,type,switch,maxInt)
    
        num = num - 10

    print("num=%d\n" %(num))
    if(num>=1):
        #----------------------------------cos_series-------------------------------------#
        type = 0
        print("Test de cos_series, maxInt=%d\n" %(maxInt))    
        ##----------------variable = nombre de termes-------------------##
        print("\tcos TERMES")
        file_abs = rootpath + "cos_terme.txt"
        reportWriter(file_abs,1+shift,niveauPrecisionFixe,nbIteration,nbRepete,cos_series,type,switch,maxInt)
        ##---------------variable = niveau de precision-----------------##
        print("\tcos PRECISION")
        file_abs = rootpath + "cos_precision.txt"
        reportWriter(file_abs,2+shift,nbTermesFixe,nbIteration,nbRepete,cos_series,type,switch,maxInt)
    
    print("num=%d\n" %(num))
    print("Fin de test\n")
    
def test2(nbIteration, switch):
    rootpath = "/usr/share/SageMath/2I013/data/"
    maxInt = 100
    #nbIteration = 100
    nbRepete = 5
    moyen = 0
    file_abs = ""
    nbTermesFixe = 12
    niveauPrecisionFixe = 20
    type = 0
    print("\texp_zero_zero PRECISION special")
    file_abs = rootpath + "exp_zero_zero_precision_special.txt"
    reportWriter(file_abs,3,nbTermesFixe,20,nbRepete,exp_zero_zero_series,type,switch,maxInt)


"""
def comp_inverse_un(switch,nbTermes,maxInt,precision):
    
    #C'est une fonction qui change le polynome random pour que la fonction 'inverse_un_series' puis le prendre comme un parametre, et effectuera cette fonction
    
    P = poly_random(switch,nbTermes,maxInt)
    P = P - P(0) + P.parent().one()
    return inverse_un_series(P,precision)


def comp_log_zero_un(switch,nbTermes,maxInt,precision):
    P = poly_random(switch,nbTermes,maxInt)
    P = P - P(0) + P.parent().one()
    return log_zero_un_series(P,precision)


def comp_exp_zero_zero(switch,nbTermes,maxInt,precision):
    P = poly_random(switch,nbTermes,maxInt)
    P = P - P(0) + P.parent().zero()
    return exp_zero_zero_series(P,precision)


def comp_arctan_zero_zero(switch,nbTermes,maxInt,precision):
    P = poly_random(switch,nbTermes,maxInt)
    P = P - P(0) + P.parent().zero()
    return arctan_series(P,precision)


def comp_tan_zero_zero(switch,nbTermes,maxInt,precision):
    P = poly_random(switch,nbTermes,maxInt)
    P = P - P(0) + P.parent().zero()
    return tan_series(P,precision)


def comp_sin_zero_zero(switch,nbTermes,maxInt,precision):
    P = poly_random(switch,nbTermes,maxInt)
    P = P - P(0) + P.parent().zero()
    return sin_series(P,precision)


def comp_cos_zero_zero(switch,nbTermes,maxInt,precision):
    P = poly_random(switch,nbTermes,maxInt)
    P = P - P(0) + P.parent().zero()
    return cos_series(P,precision)
"""
