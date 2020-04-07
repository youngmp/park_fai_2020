# script to find cusps in pi4/pi5 space for a fixed zeta

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import brentq

exp = np.exp

#forcevelocitycurvesF_A

#XXXXXXX
def fp_pos(U,pi3,pi4,pi5,pi6):
    
    return -(((1+pi3)*(exp(pi4)*(1-exp(pi5-pi5/(pi6*U)))-(1-exp(-pi5/(pi6*U)))*(1-pi6*U)))/((-1+exp(pi4))*(1+(1-exp(-pi5/(pi6*U)))*pi3)*(1-pi6*U)))

#XXXXXXX
def fp_neg(u,pi3,pi4,pi5,pi6):
    return -(1+pi6*u/(exp(pi4)-1))/(1-pi6*u)

#forcevelocitycurvesF_{-A}

#XXXXXXX
def fm_pos(u,pi3,pi4,pi5,pi6):
    return (1-pi6*u/(exp(pi4)-1))/(1+pi6*u)

#XXXXXXX
def fm_neg(U,pi3,pi4,pi5,pi6):
    return ((1+pi3)*(exp(pi4)*(1-exp(pi5+pi5/(pi6*U)))-(1-exp(pi5/(pi6*U)))*(1+pi6*U)))/((-1+exp(pi4))*(1+(1-exp(pi5/(pi6*U)))*pi3)*(1+pi6*U))

#DerivativeforcevelocitycurvesF_A

#XXXXXX
def dfp_pos(U,pi3,pi4,pi5,pi6):
    return -(((1+pi3)*(exp(pi4+pi5)*pi3*pi6**2*U**2+exp(pi4+(2*pi5)/(pi6*U))*(1+pi3)*pi6**2*U**2+exp(pi5/(pi6*U))*pi5*(-1+pi6*U)**2+exp(pi4+pi5+pi5/(pi6*U))*(1+pi3)*(-(pi6**2*U**2)+pi5*(-1+pi6*U))-exp(pi4+pi5/(pi6*U))*pi3*(pi6**2*U**2+pi5*(-1+pi6*U))))/((-1+exp(pi4))*(pi3-exp(pi5/(pi6*U))*(1+pi3))**2*pi6*U**2*(-1+pi6*U)**2))

def dfp_neg(U,pi3,pi4,pi5,pi6):
    return -((exp(pi4)*pi6)/((-1+exp(pi4))*(-1+pi6*U)**2))

#DerivativeforcevelocitycurvesF_{-A}
def dfm_pos(U,pi3,pi4,pi5,pi6):
    return -((exp(pi4)*pi6)/((-1+exp(pi4))*(1+pi6*U)**2))

#XXXXXXX
def dfm_neg(U,pi3,pi4,pi5,pi6):
    return -(((1+pi3)*(exp(pi4+pi5+(2*pi5)/(pi6*U))*pi3*pi6**2*U**2+exp(pi4)*(1+pi3)*pi6**2*U**2+exp(pi5/(pi6*U))*pi5*(1+pi6*U)**2+exp(pi4+pi5/(pi6*U))*pi3*(pi5+pi5*pi6*U-pi6**2*U**2)-exp(pi4+pi5+pi5/(pi6*U))*(1+pi3)*(pi5+pi5*pi6*U+pi6**2*U**2)))/((-1+exp(pi4))*(-1+(-1+exp(pi5/(pi6*U)))*pi3)**2*pi6*U**2*(1+pi6*U)**2))

#2ndderivativeforcevelocity

#XXXXXXX
def d2fp_pos(U,pi3,pi4,pi5,pi6):
    aa=2*exp((2*pi5)/(pi6*U))*(1+pi3)**2*pi6**4*U**4+exp(pi5/(pi6*U))*pi3*(1+pi3)*(-4*pi6**4*U**4+pi5**2*(-1+pi6*U)**2-2*pi5*pi6*U*(1-3*pi6*U+2*pi6**2*U**2))+pi3**2*(2*pi6**4*U**4+pi5**2*(-1+pi6*U)**2+2*pi5*pi6*U*(1-3*pi6*U+2*pi6**2*U**2))
    return -(((1+pi3)*(-((pi3-exp(pi5/(pi6*U))*(1+pi3))**2*pi5*(-1+pi6*U)**2*(-2*(-1+exp(pi4+pi5))*pi6*U+pi5*(-1+exp(pi4+pi5)+pi6*U)))-2*(-pi3+exp(pi5/(pi6*U))*(1+pi3))*(-1+pi6*U)*((-1+exp(pi5/(pi6*U)))*pi6**2*U**2-pi5*(-1+exp(pi4+pi5)+pi6*U))*(exp(pi5/(pi6*U))*pi6**2*U**2+pi3*(pi5-pi5*pi6*U+(-1+exp(pi5/(pi6*U)))*pi6**2*U**2))-(-1+exp(pi4+pi5)-exp(pi4+pi5/(pi6*U))+pi6*U+exp(pi5/(pi6*U))*(1-pi6*U))*(aa)))/(exp((3*pi5)/(pi6*U))*(-1+exp(pi4))*(1+pi3-pi3/exp(pi5/(pi6*U)))**3*pi6**2*U**4*(1-pi6*U)**3))

def d2fp_neg(U,pi3,pi4,pi5,pi6):
    return (2*exp(pi4)*pi6**2)/((-1+exp(pi4))*(-1+pi6*U)**3)

def d2fm_pos(U,pi3,pi4,pi5,pi6):
    return (2*exp(pi4)*pi6**2)/((-1+exp(pi4))*(1+pi6*U)**3)

#XXXXXX
def d2fm_neg(U,pi3,pi4,pi5,pi6):
    return ((1+pi3)*(exp(pi5/(pi6*U))*pi3*pi5*(1+pi6*U)**2*((1+pi3+exp(pi5/(pi6*U))*pi3)*pi5-2*(-1+(-1+exp(pi5/(pi6*U)))*pi3)*pi6*U)*(-1+exp(pi4)-exp(pi4+pi5+pi5/(pi6*U))-pi6*U+exp(pi5/(pi6*U))*(1+pi6*U))-2*exp(pi5/(pi6*U))*pi3*(-1+(-1+exp(pi5/(pi6*U)))*pi3)*pi5*(1+pi6*U)*(exp(pi4)*pi6**2*U**2+exp(pi5/(pi6*U))*pi5*(1+pi6*U)**2-exp(pi4+pi5+pi5/(pi6*U))*(pi5+pi5*pi6*U+pi6**2*U**2))+(1+pi3-exp(pi5/(pi6*U))*pi3)**2*(2*exp(pi4)*pi6**4*U**4+exp(pi5/(pi6*U))*pi5*(1+pi6*U)**3*(pi5+2*pi6*U)-exp(pi4+pi5+pi5/(pi6*U))*(2*pi6**4*U**4+(pi5+pi5*pi6*U)**2+2*pi5*pi6*U*(1+3*pi6*U+2*pi6**2*U**2)))))/((-1+exp(pi4))*(1+pi3-exp(pi5/(pi6*U))*pi3)**3*pi6**2*U**4*(1+pi6*U)**3)

"""
ignore scalar case for now
"""

def fp(u,pi3,pi4,pi5,pi6):
    out = np.zeros(len(u))

    out[u<0]=fp_neg(u[u<0],pi3,pi4,pi5,pi6)
    out[u>=0]=fp_pos(u[u>=0],pi3,pi4,pi5,pi6)

    return out

def fm(u,pi3,pi4,pi5,pi6):
    out = np.zeros(len(u))

    out[u<0]=fm_neg(u[u<0],pi3,pi4,pi5,pi6)
    out[u>=0]=fm_pos(u[u>=0],pi3,pi4,pi5,pi6)

    return out


def fnet(u,pi3,pi4,pi5,pi6):
    out = np.zeros(len(u))    
    out[u!=0] = p1*fm(u[u!=0],pi3,pi4,pi5,pi6)+(1-p1)*fp(u[u!=0],pi3,pi4,pi5,pi6)
    out[u==0] = -1+2*p1
    return out


def dfp(u,pi3,pi4,pi5,pi6):
    out = np.zeros(len(u))
    out[u<0] = dfp_neg(u[u<0],pi3,pi4,pi5,pi6)
    out[u>=0] = dfp_pos(u[u>=0],pi3,pi4,pi5,pi6)
    return out

def dfm(u,pi3,pi4,pi5,pi6):
    out = np.zeros(len(u))
    out[u<0] = dfm_neg(u[u<0],pi3,pi4,pi5,pi6)
    out[u>=0] = dfm_pos(u[u>=0],pi3,pi4,pi5,pi6)
    return out

def dfnet(u,pi3,pi4,pi5,pi6):
    out = np.zeros(len(u))
    out[u!=0]=p1*dfm(u,pi3,pi4,pi5,pi6)+(1-p1)*dfp(u,pi3,pi4,pi5,pi6)
    out[u==0]=(exp(pi4)*pi6)/(1 - exp(pi4))
    return out


def d2fp(u,pi3,pi4,pi5,pi6):
    out = np.zeros(len(u))
    out[u<0] = d2fp_neg(u[u<0],pi3,pi4,pi5,pi6)
    out[u>=0] = d2fp_pos(u[u>=0],pi3,pi4,pi5,pi6)

    return out

def d2fm(u,pi3,pi4,pi5,pi6):
    out = np.zeros(len(u))
    out[u<0] = d2fm_neg(u[u<0],pi3,pi4,pi5,pi6)
    out[u>=0] = d2fm_pos(u[u>=0],pi3,pi4,pi5,pi6)
    return out

def d2fnet(u,pi3,pi4,pi5,pi6):
    return p1*d2fm(u,pi3,pi4,pi5,pi6)+(1-p1)*d2fp(u,pi3,pi4,pi5,pi6)


# nullclines
def n1(u,z,pi3,pi4,pi5,pi6):
    return (fp(u,pi3,pi4,pi5,pi6)-z*u)/(fp(u,pi3,pi4,pi5,pi6)-fm(u,pi3,pi4,pi5,pi6))
    #a=1/(fp(u,pi3,pi4,pi5,pi6)-fm(u,pi3,pi4,pi5,pi6))
    #return exp(-np.abs(a))*a

def n2(u,z,pi3,pi4,pi5,pi6):
    return (dfp(u,pi3,pi4,pi5,pi6)-z)/(dfp(u,pi3,pi4,pi5,pi6)-dfm(u,pi3,pi4,pi5,pi6))
    #a=1/(dfp(u,pi3,pi4,pi5,pi6)-dfm(u,pi3,pi4,pi5,pi6))
    #return exp(-np.abs(a))*a



def zero_crossings(pi5,u,z,pi3,pi4,pi6):
    """
    return 1 if 2, -1 if 1 or 0
    """

    ndiff = n1(u,z,pi3,pi4,pi5,pi6)-n2(u,z,pi3,pi4,pi5,pi6)
    
    zero_crossings = np.where(np.diff(np.signbit(ndiff)))[0]#np.where(np.diff(np.sign(ndiff)))[0]

    # only take u where |ndiff|<=2
    root_u = u[zero_crossings][np.abs(ndiff[zero_crossings])<=2]
    root_n = ndiff[zero_crossings][np.abs(ndiff[zero_crossings])<=2]

    #print(pi5)
    if len(root_u) == 2:
        return 1

    else:
        return -1

def main():

    pi3=1
    pi4=4.7
    pi5=.1
    pi6=10
    z=1

    


    p56=pi5/pi6
    ep4=exp(pi4)
    ep5=exp(pi5)



    u = np.linspace(-.1,-.00,100000)

    # get crossing point(s)

    pi4_array = np.linspace(0,5,10)
    pi5_array = np.zeros(len(pi4_array))

    for i in range(len(pi4_array)):
        print(pi4_array[i])
        try:
            pi5_array[i] = brentq(zero_crossings,.01,2,args=(u,z,pi3,pi4_array[i],pi6),rtol=1e-4)
        except ValueError:
            pi5_array[i] = np.nan

    fig = plt.figure()
    ax = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    
    ax.plot(pi4_array,pi5_array)


    ndiff = n1(u,0,pi3,1.6667,.4,10)-n2(u,0,pi3,1.6667,.4,10)
    ax2.plot(ndiff)

    ax2.set_ylim(-1,1)

    plt.show()


if __name__ == '__main__':
    main()
