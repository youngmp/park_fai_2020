"""
creative space to generate movies of the lubrication model

avconv -r 40 -start_number 1 -i test%d.png -b:v 1000k test.mp4
"""

import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] 

from lubrication import lubrication as LB


pi=np.pi
exp=np.exp



def main():

    # name of movie
    movdir = 'mov'
    if not(os.path.exists(movdir)):
        os.makedirs(movdir)

    #a = LB(T=200,Z0=-5,U0=.05,phi1=0.57,pi5=0.02,Rp=1.5,Rc=2.15)
    #a = LB(T=200,Z0=-5,U0=.05,phi1=0.57,pi5=0.02,Rp=1.5,Rc=2.15)

    #a = LB(T=200,Z0=-5,U0=.05,phi1=0.5,pi5=0.1,Rp=0.96,Rc=1.22)
    a = LB(T=275,Z0=-5,U0=.05,phi1=0.5,pi5=0.02,Rp=1.5,Rc=2.15,F0=200)
    
    #a = LB(T=275,Z0=-5,U0=.05,phi1=0.5,pi5=0.02,Rp=1.5,Rc=2.15,eta_fast=3e-2)
    
    a.run_sim()

    skipn = 50

    if True:
        fig = plt.figure()
        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312)

        ax1.set_title('movie preview')
        ax1.plot(a.t,a.Z*a.Rc)
        ax2.plot(a.t,a.U*a.U_scale)

        plt.show()
        plt.close()

    # choose frames to use
    frames = np.arange(0,len(a.U),skipn)
    
    for i in range(len(frames)):
        fig = plt.figure(figsize=(5,6))
        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312)
        ax3 = fig.add_subplot(313)

        j = frames[i]
        ax1.plot(a.t[:j],a.Z[:j]*a.Rc)
        ax2.plot(a.t[:j],a.U[:j]*a.U_scale/10)

        x = np.linspace(-10*a.Rc,10*a.Rc,100)
        x2 = np.linspace(0,2*pi,20)
        
        ax3.plot(a.Rp*np.cos(x2)+a.Z[j]*a.Rc,a.Rp*np.sin(x2))
        ax3.plot(x,a.pi1(x),c='k')
        ax3.plot(x,-a.pi1(x),c='k')

        ax1.set_title(r'$\phi_1 = \frac{\text{Forward Motors}}{\text{Backward Motors}}=$'+str(a.phi1)+r'\\Vesicle Radius='+str(a.Rp)+r'$\mu$m'+r'\\Channel Radius='+str(a.Rc)+r'$\mu$m')        
        ax3.set_title('Viscous Drag=%.2e'% (a.viscous_drag(a.Z[j])*a.zeta_inf*1e-8)+'Kg/s')

        ax2.set_xlabel('Time (s)')
        ax3.set_xlabel('Position ($\mu$m)')
        
        ax1.set_ylabel('Position ($\mu$m)')
        ax2.set_ylabel('Velocity ($\mu$m/s)')
        
        ax1.set_xlim(0,a.T)
        ax2.set_xlim(0,a.T)
        ax3.set_xlim(-6*a.Rc,6*a.Rc)
        
        ax1.set_ylim(np.amin(a.Z*a.Rc)-.5,np.amax(a.Z*a.Rc)+.5)
        ax2.set_ylim(np.amin(a.U*a.U_scale)-.1,np.amax(a.U*a.U_scale)+.1)

        plt.tight_layout()

        plt.savefig(movdir+'/'+str(i)+'.png')
        plt.close()


if __name__ == "__main__":
    main()
