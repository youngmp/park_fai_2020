"""
simulation of langevin equation derived to understand Allard 2019

\begin{equation*}
 dx = -A(x)dt + B(x) dW(t),
\end{equation*}
where $A(x) = \kappa^-(x) - \kappa^+(x)$, $B(x) = \kappa^+(x)+\kappa^-(x)$, and
\begin{equation}\label{eq:continuous_rates}
\kappa^-(x) = \frac{x}{L} \kappa_0 \exp\left(\gamma \left(1-\frac{x}{L}\right)\right), \quad\quad \kappa^+(x) = \left(1-\frac{x}{L}\right)\kappa_0 \exp\left(\gamma \frac{x}{L}\right).
\end{equation}

surely many good general reference exist, but here we use Equation 10.2 from Ermentrout and Terman to numerically integrate the Langevin equation.

"""

import argparse

import os
import time
import numpy as np

import matplotlib.pyplot as plt
import scipy.integrate as integrate

from scipy.interpolate import interp1d
from mpl_toolkits.mplot3d import Axes3D

import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] 

pi=np.pi
exp=np.exp
sqrt=np.sqrt

# parameters from Allard
kap0 = 1/2
gm = 3
K = 35
#K = 15
ze = 0



def kn(x):
    return x*kap0*exp(gm*(1-x))

def kp(x):
    return (1-x)*kap0*exp(gm*x)

def kn_drag(x,ze=0,k0=.5):
    aa=1/(1+2*ze/K)
    bb=ze/(K+2*ze)
    return x*k0*exp(gm*(aa*(1-x)+bb))

def kp_drag(x,ze=0,k0=.5):
    aa=1/(1+2*ze/K)
    bb=ze/(K+2*ze)
    return (1-x)*k0*exp(gm*(aa*x+bb))

def kn_discrete(i):
    if i == 0:
        return 0
    else:
        return i*kap0*exp(gm*(1-i/K))

def kp_discrete(i):
    if i == K:
        return 0
    else:
        return (K-i)*kap0*exp(gm*(i/K))

def A(x):
    return (kp_drag(x)-kn_drag(x))/K

def B(x):
    return (kp_drag(x)+kn_drag(x))/K**2

def C(x):
    """
    noise coeffieicnt
    """
    return sqrt(B(x))

def check_boundary(x_current,x_previous,choice='reflect1'):

    if choice == 'reflect1':
        if x_current > 1:
            return 1 - (x_current - 1)

        elif x_current < 0:
            return -x_current

        else:
            return x_current
        
    elif choice == 'reflect2':
        if (x_current > 1) or (x_current < 0):
            return x_previous

    elif choice == 'sticky':        
        if x_current > 1:
            return 1

        elif x_current < 0:
            return 0

        else:
            return x_current

    else:
        return x_current


def run_langevin(t,switch_motor=None):

    dt = t[-1]-t[-2]

    x = np.zeros(len(t))
    x[0] = 0 # explicit for clarity
    
    np.random.seed(0)

    
    switch_times = []
    side = 0
    
    for i in range(len(t)-1):
        
        # Euler step
        w = np.random.normal(0,1)
        x[i+1] = x[i] + dt*( A(x[i]) ) + sqrt(dt)*C(x[i])*w

        # check boundary condition
        x[i+1] = check_boundary(x[i+1],x[i])

        # track switch times
        if switch_motor != None:
                        
            if (side == 0) and (x[i+1] >= 1-switch_motor):
                side = 1
                switch_times.append(t[i])

            if (side == 1) and (x[i+1] <= switch_motor):
                side = 0
                switch_times.append(t[i])
    
    return x,switch_times

def run_gillespie(T,i0=0,seed=0):
    t = 0
    i = i0

    #sol =
    sol = []
    t_array = []
    switch_times = []

    sol.append(i0)
    t_array.append(0)

    switch_motor = 0 # other motor is (self.N-2)
    side = 0

    np.random.seed(seed)
    while t < T:

        r1,r2 = np.random.uniform(size=(2,))
        a1 = kp_discrete(i) # propensity of growth
        a2 = kn_discrete(i) # propensity of decay
        a0 = a1+a2

        # next "reaction"
        tau = 1/a0 * np.log(1/r1)

        # increment up
        if (r2 >= 0) and (r2 <= a1/a0):
            i += 1
        if (r2 > a1/a0) and (r2 <= 1):
            i -= 1

        t += tau
        
        
        if (side == 0) and (i >= K-switch_motor):
            side = 1
            switch_times.append(t)

        if (side == 1) and (i <= switch_motor):
            side = 0
            switch_times.append(t)


        sol.append(i)
        t_array.append(t)
    
        
    return t_array,sol,switch_times,t

def fp_ss(x):
    """
    fokker-planck steady-state

    x: must be array
    """
    dx = x[-1]-x[-2]

    integral = integrate.cumtrapz( (A(x)/B(x)),x, initial=0 )
    
    # normalization factor
    c = np.sum(np.exp(2*integral)/B(x))*dx
    
    return np.exp(2*integral)/(c*B(x))

def psi(x):
    """
    function psi, used in probability/first passage time problems, defined in Gardiner, for convenience.

    x: domain vector must be array
    zeta: scalar
    """

    #dx = (domain[-1]-domain[0])/len(domain)

    integral = integrate.cumtrapz( (A(x)/B(x)),x, initial=0 )
    
    f = exp( 2*integral )

    return f
    #return interp1d(x,f)(x)


def pi_a(x,d=False):
    """
    function pi_a, probability of escale through (lower) boundary a
    x must be array
    """

    if d:
        return -pi_b(x,domain,d=True)
    else:
        return 1 - pi_b(x,domain)


def pi_b(x,d=False):

    dx = (domain[-1]-domain[-2])/len(domain)
    tot = np.sum(psi(domain,domain))*dx
    f = np.cumsum(psi(domain,domain))*dx

    if d:
        return psi(x,domain)/tot
    
    else:
        return interp1d(domain,f)(x)/tot

def mfpt(x,domain):
    """
    in the case of one reflecting boundary and one absorbing boundary (as in the case of allard), we use formula 5.2.160 in Gardiner.
    """
    
    a = 0 # reflecting boundary is on left side
    dx = (domain[-1]-domain[0])/len(domain)

    p = psi(domain)
    f1 = integrate.cumtrapz( p/B(domain),domain, initial=0 )

    integral = integrate.cumtrapz( f1/p,domain, initial=0 )

    if False:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        #ax.plot(domain,psi(domain,domain)/B(domain)*dx)
        ax.plot(domain,f1)
        plt.show()
    
    tot = np.sum(f1/p)*dx
    
    return 2*(tot-integral)
    
def mfpt_ode(y,t,domain=np.zeros(1),zeta=0):
    """
    for the calculation of mean first passage time
    2nd order ode
    this formula is useful when there are two absorbing boundaries as in the case of Myosin mean-field.
    """

    y1 = y[0]
    y2 = y[1]

    dx = (domain[-1]-domain[-2])/len(domain)
    tot = np.sum(psi(domain,domain))*dx

    A_val = A(t)
    pi_a_val = pi_a(t,domain=domain)
    dpi_a_val = pi_a(t,domain=domain,d=True)

    B_val = B(t)

    return np.array([y2,2*(-pi_a_val-2*y1-(A_val*pi_a_val + B_val*dpi_a_val)*y2)/(B_val*dpi_a_val)])


def run_mfpt(x_test):

    # pick a zeta parameter. pi_a,b, psi will be computed once for a given zeta, then recomputed (without history) for another zeta.
    
    zeta_par = 1
    #x_test = np.linspace(-9e-2,9e-2,201)
    dx = (x_test[-1]-x_test[0])/len(x_test)

    mfpt = np.zeros((len(x_test),2))
    zeta = 1

    for i in range(len(x_test)-1):
        t = x_test[i]
        mfpt[i+1,:] = mfpt[i,:] + dx*mfpt_ode(mfpt[i,:],t,domain=x_test,zeta=zeta_par)

    return mfpt

def plot_mfpt(peak_pos):
    """
    make sure inputs to this plot have 2 peaks
    """
    print('peak_pos',peak_pos)

    fig = plt.figure()
    fig.set_figheight(3)
    fig.set_figwidth(8)

    ax1 = fig.add_subplot(141)
    ax2 = fig.add_subplot(142)
    ax3 = fig.add_subplot(143)
    ax4 = fig.add_subplot(144)


    x_test = np.linspace(peak_pos[0],peak_pos[1],201)
    #mfpt = run_mfpt(x_test)
    mfpt_val = mfpt(x_test,x_test)
    
    #print('mfpt_val',mfpt_val)
    
    ax1.plot(x_test,psi(x_test))
    #ax2.plot(x_test,pi_a(x_test,x_test))
    #ax3.plot(x_test,pi_b(x_test,x_test))
    ax4.plot(x_test,mfpt_val)

    ax1.set_title(r'$\psi(x)$')
    #ax2.set_title(r'$\pi_a(x)$')
    #ax3.set_title(r'$\pi_b(x)$')

    ax2.set_xlabel(r'$x$')
    ax2.set_xlabel(r'$x$')

    ax2.set_ylabel('Probably of Escape (Left)')
    ax3.set_ylabel('Probably of Escape (Right)')
    ax4.set_ylabel('MFPT (Right)')


    plt.tight_layout()


def main():

    parser = argparse.ArgumentParser(description='run the agent-based Myosin motor model',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-s','--seed',default=0,type=int,help='Set random number seed for simulation')
    #parser.add_argument('-z','--zeta',default=0,type=np.float64,help='Set viscous drag')

    #parser.add_argument('-N','--N',default=35,type=int,help='Set total motor number')
    #parser.add_argument('-Y','--nY',default=100,type=int,help='Set total motor number (Y right preferrred)')

    parser.add_argument('--mfpt_fname',default='',type=str,help='If str is defined, then the simulation will save mfpt data to the directory')

    args = parser.parse_args()
    print(args)
    
    Tfinal = 30000
    t_gillespie,sol_gillespie,switch_times,total_time = run_gillespie(Tfinal,seed=args.seed)

    list_of_durations = np.diff(switch_times)
    list_of_durations = np.append(np.array([total_time]),list_of_durations)

    if args.mfpt_fname != '':
        np.savetxt(args.mfpt_fname,list_of_durations)

    # get peak positions    
    dom = np.linspace(0,1,101)
    a = fp_ss(dom)
    peak_idxs = np.r_[True, a[1:] > a[:-1]] & np.r_[a[:-1] > a[1:], True]
    
    #print(peak_idxs.dtype)
    peak_pos = dom[peak_idxs]
    plot_mfpt(peak_pos)

    if True:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        #print('A(0)/B(0)+A(1)/B(1)',A(0)/B(0)+A(1)/B(1))

        #ratio = A(dom)/B(dom)
        #ax.plot(dom,ratio)
        #print('ratio[0],ratio[-1]',ratio[0],ratio[-1])
        
        #ax.plot(dom,B(dom))

        dx = (dom[-1]-dom[0])/len(dom)#dom[-1]-dom[-2]

        integral = np.cumsum(A(dom)/B(dom))*dx
        ax.plot(dom,integral)
        ax.plot(dom,fp_ss(dom))
        #ax.plot(np.exp(2*integral)/B(x))
        #plt.show()


    

    
    fig = plt.figure()
    ax11 = fig.add_subplot(221)
    ax12 = fig.add_subplot(222)
    ax21 = fig.add_subplot(223)
    ax22 = fig.add_subplot(224)


    # run langevin
    dt = 0.01 # time step for langevin simulation
    t_langevin = np.linspace(0,Tfinal,int(Tfinal/dt))
    sol_langevin,switch_times = run_langevin(t_langevin,switch_motor=peak_pos[0])

    #print(switch_times)
    print('MFPT Langevin',np.mean(np.diff(switch_times)))
    print('Total switches in Langevin',len(np.diff(switch_times)))
    
    
    # run Gillespie algorithm
    ax11.plot(t_gillespie,sol_gillespie)

    #https://stackoverflow.com/questions/3866520/plotting-histograms-whose-bar-heights-sum-to-1-in-matplotlib/16399202#16399202
    #weights = np.ones_like(sol_gillespie)/len(sol_gillespie)
    #ax12.hist(sol_gillespie,bins=int(K*1.),density=True) # plot distribution of agent-based model
    

    dom = np.linspace(0,1,100)


    ax21.plot(t_langevin,sol_langevin) # plot langevin solution
    ax21.scatter(switch_times,np.zeros(len(switch_times)),color='tab:red')
    
    ax22.hist(sol_langevin,bins=100,density=True) # plot langevin distribution (histogram of states of one solution)

    

    
    ax22.plot(dom,fp_ss(dom),color='r') # plot FP steady state    
    ax11.set_ylabel(r'State')
    
    ax21.set_xlabel(r'$t$')
    ax21.set_ylabel(r'State')

    ax22.set_xlabel(r'State')
    ax12.set_ylabel(r'Count')
    ax22.set_ylabel(r'Count')

    ax11.set_title('Gillespie')

    ax21.set_title('Langevin')

    plt.show()


    #plt.tight_layout()



    #fig2 = plt.figure()
    #bx = fig2.add_subplot(111)

    
    #plt.show()

if __name__ == "__main__":
    main()
