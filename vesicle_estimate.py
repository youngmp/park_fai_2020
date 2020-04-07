import numpy as np

from lubrication import lubrication

def main():
    
    # initialize llubrication model with appropriate height.
    # assume constant constriction diameter
    # check yptgf1/vesicle_estimates/vesicle_figure.svg and .png as wel as vesicle_estimate1.txt

    # in Silva et al Kapitein 2015 fig 1, height is roughly 0.08um
    # in fig2, hight at constriction is rougly 0.15um, height in head is roughly 0.5um. Take Rp=1 and Rc=Rp+height.
    a = lubrication(phi1=.57,Rp=0.5,Rc=.65,
                    pi3=1,pi4=4.7,pi5=0.1,pi6=10,
                    mu=1.2,T=1,constriction='constant',U0=0.2,
                    dt=0.02,eps=1,
                    F0=50)
    
    ze_scale1 = 1e-6*(6*np.pi*a.Rp*a.mu) # kg/s
    ze = a.viscous_drag(0)*ze_scale1

    print(ze)

    
    
if __name__ == "__main__":
    main()
