# figure generation file
# all bifurcation diagrams are generated using XPPAUTO. 

from decimal import Decimal
from matplotlib.collections import PatchCollection

import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, mark_inset)


from lib import collect_disjoint_branches
from lubrication import lubrication

from bifurcations.cusp import *

import matplotlib as mpl

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['text.usetex'] = True
mpl.rcParams['pgf.texsystem'] = 'pdflatex'
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}',
                                       r'\usepackage{siunitx}']

# font size
size = 12

exp = np.exp
pi = np.pi

red = '#e66101'
yellow = '#fdb863'
lightpurp = '#b2abd2'
purp = '#5e3c99'


# color regions. 1 stable => green. 3 stable => yellow. 5 stable => blue
# same colors used in inkscape
color1pts = '#00ff00';color1pts_al = .15
color3pts = '#ffffbd';color3pts_al = .47
color5pts = '#0000ff';color5pts_al = .18



class Arrow3D(FancyArrowPatch):
    """
    A class for drawing arrows in 3d plots.
    """
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)


def fill_axis(ax_main,keys,val_dict,type_dict,
              ze_scale,tol=1e-3,beautify=True):
    """
    this function is to clean up the 2 parameter bifuraction data. most
    things here are hard-coded and will not generalize very well.
    I put all this junk in this function so it can be accessed again for
    a minimal version of this data in another figure.
    """
    for i in range(len(keys)):
        err1 = np.sum(np.abs(np.diff(val_dict[keys[i]][:,0])))
        err2 = np.sum(np.abs(np.diff(val_dict[keys[i]][:,1])))

        label = ''

        zorder = 2
        # this is a 2 par diagram
        # so ignore keys where one coordinate does not change.
        if (err1>tol) and (err2>tol) and (keys[i] != 'br22')\
            and (keys[i] != 'br34'):

            markersize = 15
            if keys[i] == 'br29':
                marker = '*'
                skipn = 200
                label = ''
                color = red
                #color = '#19439c'
                mew = 1
                markersize = 20

            elif keys[i] == 'br28':
                marker = ''
                skipn = 200
                label = ''
                color = red
                #color = '#19439c'
                mew = 1
                markersize = 20

                
            elif (keys[i] == 'br27'):
                marker = '+'
                skipn = 200
                label = ''
                color = yellow
                #color = '#535f42'
                mew = 5
                markersize = 20

            elif keys[i] == 'br26':
                marker = ''
                skipn = 200
                label = ''
                color = yellow
                #color = '#535f42'
                mew = 5
                markersize = 20

                
            elif keys[i] == 'br31':
                marker = 'x'
                skipn = 240
                label = ''
                color = lightpurp
                #color = '#997f41'
                mew = 5

            elif keys[i] == 'br30':
                marker = ''
                skipn = 240
                label = ''
                color = lightpurp
                #color = '#997f41'
                mew = 5

            elif keys[i] == 'br32':
                marker = ''
                skipn = 200
                label = ''
                color = purp
                #color = '#3f0000'
                mew = 5
                markersize = 20
                zorder = 3

            elif keys[i] == 'br33':
                marker = '4'
                skipn = 200
                label = ''
                color = purp
                #color = '#3f0000'
                mew = 5
                markersize = 20
                zorder = 3
                
            else:
                marker = ''

            kw = {'color':color,'lw':0,'marker':marker,
                  'markersize':markersize,'mew':mew,'zorder':zorder}
            
            x_vals = val_dict[keys[i]][:,0]
            y_vals = val_dict[keys[i]][:,1]*ze_scale
            
            print('valid i=',i,'errs',err1,err2,'key=',keys[i])


            if not(beautify):
                color='k'
                label=''
                kw = {'color':'k','lw':0}
                
            ax_main.plot(x_vals,y_vals,
                         color=color,label=label,lw=3,zorder=i)
                
            if marker != '':
                if keys[i] == 'br27':
                    cut = 7000
                    # need to break this up because point density
                    # increases halfway through the curve
                    ax_main.plot(x_vals[cut::skipn],y_vals[cut::skipn],**kw)

                    ax_main.plot(x_vals[:cut:skipn*10],
                                 y_vals[:cut:skipn*10],**kw)

                else:
                    ax_main.plot(x_vals[::skipn],y_vals[::skipn],**kw)

            label = ''

    return ax_main


def constriction():
    """
    plot constriction diagram with time plots
    """

    debug_time = False
    if debug_time:
        T1 = .1
        T2 = .1
    else:
        T1 = 3000
        T2 = 300
        
        
    fig = plt.figure(figsize=(8,3.3))
    gs = gridspec.GridSpec(ncols=3,nrows=2)

    ax11 = fig.add_subplot(gs[0,0])
    ax12 = fig.add_subplot(gs[0,1])
    ax13 = fig.add_subplot(gs[0,2])

    ax21 = fig.add_subplot(gs[1,0])
    ax22 = fig.add_subplot(gs[1,1])
    ax23 = fig.add_subplot(gs[1,2])

    ################################################################ example 1
    # T=3000
    a = lubrication(phi1=.57,Rp=.96,Rc=1.22,
                    pi3=1,pi4=4.7,pi5=0.1,pi6=10,
                    mu=1.2,T=T1,constriction='piecewise',U0=0.2,
                    dt=0.02,eps=1,
                    F0=50,method='euler')
    a.Z0 = -5/a.Rp
    
    # run example trajectory
    a.run_sim()
    t1 = a.t*a.tau0
    U1 = a.U*a.U_scale
    Z1 = a.Z*a.Z_scale

    # run another example
    a.Z0 = -5/a.Z_scale
    a.U0 = -.3/a.U_scale
    a.T = T1
    a.run_sim()
    t1b = a.t*a.tau0
    U1b = a.U*a.U_scale
    Z1b = a.Z*a.Z_scale
    
    x = np.linspace(-7,7,1000)  # dimensional
    x2 = np.linspace(0,2*pi-(2*pi)/200,200)

    # plot vesicle
    ax11.plot(a.Rp*np.cos(x2),a.Rp*np.sin(x2),color='k',alpha=1)
    
    ax11.plot(x,a.pi1(x),color='k')
    ax11.plot(x,-a.pi1(x),color='k')
    ax11.plot([x[-1],x[-1]],[-a.base_radius,a.base_radius],color='k')

    # annotate Rp length
    ax11.annotate(r'$2R_p$', xy=(0, -a.Rc-.5),
                  xytext=(0, -a.Rc-1.5), xycoords='data', 
                  fontsize=10, ha='center', va='top',
                  arrowprops=dict(arrowstyle=('-[, widthB='
                                              +str(a.Rp)
                                              +', lengthB=.5'), lw=.5))

    # annotate Rc length channel
    arrowstyle=('-[, widthB='+str(a.Rc)+', lengthB=.5')
    ax11.annotate(r'$2R_c$', xy=(3.5, 0.0), xytext=(4, 0), xycoords='data', 
            fontsize=10, ha='left', va='center',
            arrowprops=dict(arrowstyle=arrowstyle,lw=.5))

    # annotate base diameter
    ax11.text(-7.5,1,r'$'+str(2*a.base_radius)+r'$ \si{\um}')

    # annotate constriction length
    arrowstyle = ('-[, widthB='+str(a.inner_width/2*.9)+', lengthB=.5')
    ax11.annotate(r'$'+str(a.inner_width)+r'$ \si{\um}',
                  xy=(0, a.Rc+1), xytext=(0, a.Rc+2), xycoords='data', 
                  fontsize=10, ha='center', va='bottom',
                  arrowprops=dict(arrowstyle=arrowstyle,lw=.5))

    # label channel base and starting point              
    ax11.plot([-5,-5],[-a.base_radius,a.base_radius],ls='--',color='gray')
    
    # draw arrow of motion
    ax11.annotate('',xy=(-5,0),xytext=(-a.Rp,0),
                  xycoords='data',
                  arrowprops=dict(arrowstyle='<-',facecolor='k',lw=.5))

    ax12.plot(t1[:-1],U1[:-1],color='black',label='Initial Condition 1')
    ax13.plot(t1[:-1],Z1[:-1],color='black')

    ax12.plot(t1b[:-1],U1b[:-1],color='tab:red',label='Initial Condition 2')
    ax13.plot(t1b[:-1],Z1b[:-1],color='tab:red')

    ################################################################ example 2    
    a.__init__(U0=0.02,T=T2,dt=0.02,eps=1,
               Rp=1.5,Rc=2.15,
               phi1=0.54,F0=200,n0=200,method='euler',
               pi3=1,pi4=4.7,pi5=0.02,pi6=10,constriction='piecewise')
    
    a.Z0=-5/a.Rp

    # run sim
    a.run_sim()
    t2 = a.t*a.tau0
    U2 = a.U*a.U_scale
    Z2 = a.Z*a.Z_scale

    # run another example
    
    a.Z0 = 0
    a.U0 = -.1/a.U_scale
    a.run_sim()
    t2b = a.t*a.tau0
    U2b = a.U*a.U_scale
    Z2b = a.Z*a.Z_scale
    
    # plot vesicle
    ax21.plot(a.Rp*np.cos(x2),a.Rp*np.sin(x2),color='k')

    # plot spine wall
    ax21.plot(x,a.pi1(x),color='k')
    ax21.plot(x,-a.pi1(x),color='k')
    ax21.plot([x[-1],x[-1]],[-a.base_radius,a.base_radius],color='k')
    
    # plot time traces
    ax22.plot(t2[:-1],U2[:-1],color='black')
    ax23.plot(t2[:-1],Z2[:-1],color='black')

    ax22.plot(t2b[:-1],U2b[:-1],color='tab:red')
    ax23.plot(t2b[:-1],Z2b[:-1],color='tab:red')


    ax11.axis('off')
    ax21.axis('off')

    ax22.set_xlabel(r'$t$ (s)')
    ax23.set_xlabel(r'$t$ (s)')

    ax12.set_ylabel(r'$U$ (\si{\um/s})')
    ax13.set_ylabel(r'$Z$ (\si{\um})')
    
    ax22.set_ylabel(r'$U$ (\si{\um/s})')
    ax23.set_ylabel(r'$Z$ (\si{\um})')
    
    ax11.set_ylim(x[0]/2,x[-1]/2)
    ax21.set_ylim(x[0]/2,x[-1]/2)


    ax11.set_title(r'\textbf{A}',x=0)
    ax12.set_title(r'\textbf{B}',x=0)
    ax13.set_title(r'\textbf{C}',x=0)
    
    ax21.set_title(r'\textbf{D}',x=0)
    ax22.set_title(r'\textbf{E}',x=0)
    ax23.set_title(r'\textbf{F}',x=0)

    
    ax12.legend(prop={"size":8},frameon=False,loc='lower center')

    plt.tight_layout()
    return fig
    

def critical_manifold_with_ze():
    """
    plot examples of the critical manifold.
    
    generate contour plots of the surface F(U) - zeta*U
    extract zero contour.
    """
    
    fig = plt.figure(figsize=(8,4))
    gs = gridspec.GridSpec(ncols=3,nrows=2)
    
    ax11 = fig.add_subplot(gs[0,0])
    ax12 = fig.add_subplot(gs[0,1])
    ax13 = fig.add_subplot(gs[0,2])

    ax21 = fig.add_subplot(gs[1,0])
    ax22 = fig.add_subplot(gs[1,1])
    ax23 = fig.add_subplot(gs[1,2])

    a = lubrication(phi1=.57,Rp=.96,Rc=1.22,
                    pi3=1,pi4=4.7,pi5=0.1,pi6=10,
                    mu=1.2,T=3,constriction='piecewise',U0=0.2,
                    dt=0.02,eps=1,
                    F0=50)

    ze_scale1 = 1e-6*(6*np.pi*a.Rp*a.mu) # kg/s
    
    M = 200 # total number of grid points in mesh
    
    ##### critical manifold
    U = np.linspace(-.4,.4,M)/a.U_scale # velocity range
    Z = np.linspace(-5,5,M)/a.Z_scale # position range
    Zm,Um = np.meshgrid(Z,U)
    ze = np.zeros((M,M))

    for i in range(M):
        for j in range(M):
            ze[i,j] = a.viscous_drag(Zm[i,j])
    
    C_0 = a.F(Um)-ze*Um # critical manifold

    ##### bifurcation diagram
    ze_range = np.linspace(0,20,M) # zeta range
    ze_rangem,Um = np.meshgrid(ze_range,U)
    D = a.F(Um)-ze_rangem*Um
    
    ##### plot of zeta as a function of position
    ze_plot = np.zeros(len(Z))
    
    for i in range(M):
        ze_plot[i] = a.viscous_drag(Z[i])


    # semi automating the labeling of stable/unstable curves
    ax11_contour = ax11.contour(Zm*a.Z_scale,Um*a.U_scale,C_0,[0],colors='k')
    ax11.clear()

    # ax11_contour.collections[0].get_paths()[0]
    # 0 is right curve, 1 is left, 2 is upper.
    # need to separate by min/max y value in 0 and 1 respectively.
    p0=ax11_contour.collections[0].get_paths()[0]
    p1=ax11_contour.collections[0].get_paths()[1]
    p2=ax11_contour.collections[0].get_paths()[2]
    
    x0=p0.vertices[:,0];y0=p0.vertices[:,1]
    x1=p1.vertices[:,0];y1=p1.vertices[:,1]
    x2=p2.vertices[:,0];y2=p2.vertices[:,1]

    # find x min,max of p0,p1 and use index to mark y value
    idx0 = np.argmin(x0)
    idx1 = np.argmax(x1)

    x_right_s,y_right_s = x0[:idx0],y0[:idx0] # right stable branch
    x_right_u,y_right_u = x0[idx0:],y0[idx0:] # right unstable branch    

    x_left_s,y_left_s = x1[idx1:],y1[idx1:] # left stable branch
    x_left_u,y_left_u = x1[:idx1],y1[:idx1] # left unstable branch

    ax11.plot(x_right_s,y_right_s,color='k',
              label='Stable Manifold',lw=2)
    ax11.plot(x_right_u,y_right_u,color='tab:red',
              label='Unstable Manifold',lw=2)

    ax11.plot(x_left_s,y_left_s,color='k',lw=2)
    ax11.plot(x_left_u,y_left_u,color='tab:red',lw=2)

    ax11.plot(x2,y2,color='k',lw=2)

    ###################### plot direction of motion on stable manifolds
    head_width=0.3; head_length=0.6
    head_width_fast=0.2; head_length_fast=0.4

    
    # plot on left stable manifold
    f = 1/2;n = len(x_left_s)
    arrowstyle=('-|>,head_width='+str(head_width)
                +',head_length='+str(head_length))
    arrowdict = dict(arrowstyle=arrowstyle,color='k')
    ax11.annotate('',
                  xy=(x_left_s[int(n*f)+1],y_left_s[int(n*f)+1]),
                  xytext=(x_left_s[int(n*f)],y_left_s[int(n*f)]),
                  arrowprops=arrowdict)

    # plot on right stable manifold
    f = 1/2;n = len(x_right_s)
    ax11.annotate('',
                  xy=(x_right_s[int(n*f)+1],y_right_s[int(n*f)+1]),
                  xytext=(x_right_s[int(n*f)],y_right_s[int(n*f)]),
                  arrowprops=arrowdict)
    
    # plot on upper stable manifold
    n = len(x2)

    f = 1/10
    ax11.annotate('',
                  xy=(x2[int(n*f)],y2[int(n*f)]),
                  xytext=(x2[int(n*f)+1],y2[int(n*f)+1]),
                  arrowprops=arrowdict)

    f = 4.8/10
    ax11.annotate('',
                  xy=(x2[int(n*f)],y2[int(n*f)]),
                  xytext=(x2[int(n*f)+1],y2[int(n*f)+1]),
                  arrowprops=arrowdict)
    
    f = 8/10
    ax11.annotate('',
                  xy=(x2[int(n*f)],y2[int(n*f)]),
                  xytext=(x2[int(n*f)+1],y2[int(n*f)+1]),
                  arrowprops=arrowdict)

    ####################### plot gray lines showing fast dynamics    
    # plot arrows pointing away from left stable manifold
    f=.25;n=len(x_left_u)
    arrowstyle = ('->,head_width='+str(head_width_fast)
                  +',head_length='+str(head_length_fast))
    
    arrowdict = dict(arrowstyle=arrowstyle,color='gray',ls='--')
    
    ax11.annotate('',
                  xy=(x_left_u[int(n*f)],y_left_u[int(n*f)]+.4),
                  xytext=(x_left_u[int(n*f)],y_left_u[int(n*f)]),
                  arrowprops=arrowdict)
    
    ax11.annotate('',
                  xy=(x_left_u[int(n*f)],y_left_u[int(n*f)]-.2),
                  xytext=(x_left_u[int(n*f)],y_left_u[int(n*f)]),
                  arrowprops=arrowdict)

    f=.8;
    ax11.annotate('',
                  xy=(x_left_u[int(n*f)],y_left_u[int(n*f)]+.25),
                  xytext=(x_left_u[int(n*f)],y_left_u[int(n*f)]),
                  arrowprops=arrowdict)

    ax11.annotate('',
                  xy=(x_left_u[int(n*f)],y_left_u[int(n*f)]-.1),
                  xytext=(x_left_u[int(n*f)],-.4),
                  arrowprops=arrowdict)

    
    # plot arrows pointing towards upper stable manifold
    f=.7;n = len(x2)
    ax11.annotate('',
                  xy=(x2[int(n*f)],y2[int(n*f)]+.05),
                  xytext=(x2[int(n*f)],.4),
                  arrowprops=arrowdict)


    ax11.annotate('',
                  xy=(x2[int(n*f)],y2[int(n*f)]-.05),
                  xytext=(x2[int(n*f)],-.4),
                  arrowprops=arrowdict)

    # plot arrows pointing towards upper stable manifold
    f=.5;n = len(x2)
    ax11.annotate('',
                  xy=(x2[int(n*f)],y2[int(n*f)]+.01),
                  xytext=(x2[int(n*f)],.4),
                  arrowprops=arrowdict)


    ax11.annotate('',
                  xy=(x2[int(n*f)],y2[int(n*f)]-.01),
                  xytext=(x2[int(n*f)],-.4),
                  arrowprops=arrowdict)



    # do the same as above gray arrows but for the other side
    # plot arrows pointing away from left stable manifold
    f=.75;n=len(x_right_u)
    ax11.annotate('',
                  xy=(x_right_u[int(n*f)],y_right_u[int(n*f)]+.4),
                  xytext=(x_right_u[int(n*f)],y_right_u[int(n*f)]),
                  arrowprops=arrowdict)
    
    ax11.annotate('',
                  xy=(x_right_u[int(n*f)],y_right_u[int(n*f)]-.2),
                  xytext=(x_right_u[int(n*f)],y_right_u[int(n*f)]),
                  arrowprops=arrowdict)

    f=.2;
    ax11.annotate('',
                  xy=(x_right_u[int(n*f)],y_right_u[int(n*f)]+.25),
                  xytext=(x_right_u[int(n*f)],y_right_u[int(n*f)]),
                  arrowprops=arrowdict)

    ax11.annotate('',
                  xy=(x_right_u[int(n*f)],y_right_u[int(n*f)]-.1),
                  xytext=(x_right_u[int(n*f)],-.4),
                  arrowprops=arrowdict)

    
    # plot arrows pointing towards upper stable manifold
    f=.3;n = len(x2)
    ax11.annotate('',
                  xy=(x2[int(n*f)],y2[int(n*f)]+.05),
                  xytext=(x2[int(n*f)],.4),
                  arrowprops=arrowdict)


    ax11.annotate('',
                  xy=(x2[int(n*f)],y2[int(n*f)]-.05),
                  xytext=(x2[int(n*f)],-.4),
                  arrowprops=arrowdict)

    # plot arrows pointing away from left stable manifold
    f=.25;n=len(x_left_u)
    ax11.annotate('',
                  xy=(x_left_u[int(n*f)],y_left_u[int(n*f)]+.4),
                  xytext=(x_left_u[int(n*f)],y_left_u[int(n*f)]),
                  arrowprops=arrowdict)
    
    ax11.annotate('',
                  xy=(x_left_u[int(n*f)],y_left_u[int(n*f)]-.2),
                  xytext=(x_left_u[int(n*f)],y_left_u[int(n*f)]),
                  arrowprops=arrowdict)


    # repeat stable/unstable labeling for bifurcation curve
    ax12_contour = ax12.contour(ze_rangem*ze_scale1,
                                Um*a.U_scale,D,[0],colors='k')
    ax12.clear()

    p0 = ax12_contour.collections[0].get_paths()[0]
    p1 = ax12_contour.collections[0].get_paths()[1]
    
    x0 = p0.vertices[:,0]
    y0 = p0.vertices[:,1]
    
    x1 = p1.vertices[:,0]
    y1 = p1.vertices[:,1]


    idx0 = np.argmax(x0)
    
    ax12.plot(x0[:idx0],y0[:idx0],color='tab:red',lw=2)
    ax12.plot(x0[idx0:],y0[idx0:],color='k',lw=2)
    
    ax12.plot(x1,y1,color='k',lw=2)

    #ax12.xaxis.major.formatter._useMathText = True
    
    ax13.plot(Z*a.Z_scale,ze_plot*ze_scale1,color='black',lw=2)
    
    #a = lubrication(phi1=.57,Rp=.96,Rc=1.22,pi5=.1,pi6=10,
    #                F0=50,n0=200,constriction='piecewise')
    # change parameters
    
    a.__init__(U0=0.02,T=1,dt=0.02,eps=1,
               Rp=1.5,Rc=2.15,
               phi1=0.54,F0=200,n0=200,
               pi3=1,pi4=4.7,pi5=0.02,pi6=10,constriction='piecewise')

    ze_scale2 = 1e-6*(6*np.pi*a.Rp*a.mu) # kg/s
    #ze_scale2 = 6*np.pi*1e-6 # kg/s
    #a.Rp=1.5;a.Rc=2.15;a.pi5=0.02;a.F0=200
    
    ##### critical manifold
    #U = np.linspace(-.2,.2,M) # velocity range
    #Z = np.linspace(-5,5,M) # position range
    #Zm,Um = np.meshgrid(Z,U)
    U = np.linspace(-.4,.4,M)/a.U_scale # velocity range
    Z = np.linspace(-5,5,M)/a.Z_scale # position range
    Zm,Um = np.meshgrid(Z,U)

    ze2 = np.zeros((M,M))

    for i in range(M):
        for j in range(M):
            ze2[i,j] = a.viscous_drag(Zm[i,j])
    
    C_02 = a.F(Um)-ze2*Um # critical manifold

    ##### bifurcation diagram
    #ze_range = np.linspace(0,20,M) # zeta range
    #ze_rangem,Um = np.meshgrid(ze_range,U)
    D2 = a.F(Um)-ze_rangem*Um
    
    ##### plot of zeta as a function of position
    ze_plot2 = np.zeros(len(Z))
    
    for i in range(M):
        ze_plot2[i] = a.viscous_drag(Z[i])
    

    ax21_contour = ax21.contour(Zm*a.Z_scale,Um*a.U_scale,C_02,[0],colors='k')
    ax21.clear()
    # ax21_contour.collections[0].get_paths()[0]
    # 0 is right curve, 1 is left, 2 is upper. need to separate by
    # min/max y value in 0 and 1 respectively.
    p0=ax21_contour.collections[0].get_paths()[0]
    p1=ax21_contour.collections[0].get_paths()[1]
    p2=ax21_contour.collections[0].get_paths()[2]
    
    x0=p0.vertices[:,0];y0=p0.vertices[:,1]
    x1=p1.vertices[:,0];y1=p1.vertices[:,1]
    x2=p2.vertices[:,0];y2=p2.vertices[:,1]
    
    ax21.plot(x0,y0,color='k',lw=2) # lower 
    ax21.plot(x1,y1,color='tab:red',lw=2) # middle
    ax21.plot(x2,y2,color='k',lw=2) # upper


    ####################### plot arrows showing slow dynamics
    # plot on lower
    f = .9/4;n = len(x0)
    arrowstyle = ('-|>,head_width='+str(head_width)
                  +',head_length='+str(head_length))
    arrowdict = dict(arrowstyle=arrowstyle,color='k')
    
    ax21.annotate('',xy=(x0[int(n*f)+1],y0[int(n*f)+1]),
                  xytext=(x0[int(n*f)],y0[int(n*f)]),
                  arrowprops=arrowdict)
    f = 2.1/4
    ax21.annotate('',xy=(x0[int(n*f)+1],y0[int(n*f)+1]),
                  xytext=(x0[int(n*f)],y0[int(n*f)]),
                  arrowprops=arrowdict)
    f = 3.5/4
    ax21.annotate('',xy=(x0[int(n*f)+1],y0[int(n*f)+1]),
                  xytext=(x0[int(n*f)],y0[int(n*f)]),
                  arrowprops=arrowdict)

    # plot on upper
    f = .4/4;n = len(x2)
    ax21.annotate('',xy=(x2[int(n*f)],y2[int(n*f)]),
                  xytext=(x2[int(n*f)+1],y2[int(n*f)+1]),
                  arrowprops=arrowdict)
    f = 1.9/4
    ax21.annotate('',xy=(x2[int(n*f)],y2[int(n*f)]),
                  xytext=(x2[int(n*f)+1],y2[int(n*f)+1]),
                  arrowprops=arrowdict)
    f = 3.2/4
    ax21.annotate('',xy=(x2[int(n*f)],y2[int(n*f)]),
                  xytext=(x2[int(n*f)+1],y2[int(n*f)+1]),
                  arrowprops=arrowdict)

    ####################### plot gray lines showing fast dynamics
    # plot arrows pointing away from middle unstable manifold

    f=.05;n=len(x1)
    arrowstyle = ('->,head_width='+str(head_width_fast)
                  +',head_length='+str(head_length_fast))
    arrowdict = dict(arrowstyle=arrowstyle,color='gray',ls='--')
    
    ax21.annotate('',xy=(x1[int(n*f)],y1[int(n*f)]+.4),
                  xytext=(x1[int(n*f)],y1[int(n*f)]),
                  arrowprops=arrowdict)
    
    ax21.annotate('',xy=(x1[int(n*f)],y1[int(n*f)]-.3),
                  xytext=(x1[int(n*f)],y1[int(n*f)]),
                  arrowprops=arrowdict)

    f=.25;n=len(x1)
    ax21.annotate('',xy=(x1[int(n*f)],y1[int(n*f)]+.15),
                  xytext=(x1[int(n*f)],y1[int(n*f)]),
                  arrowprops=arrowdict)
    
    ax21.annotate('',xy=(x1[int(n*f)],y1[int(n*f)]-.1),
                  xytext=(x1[int(n*f)],y1[int(n*f)]),
                  arrowprops=arrowdict)

    f=.75;n=len(x1)
    ax21.annotate('',xy=(x1[int(n*f)],y1[int(n*f)]+.15),
                  xytext=(x1[int(n*f)],y1[int(n*f)]),
                  arrowprops=arrowdict)
    
    ax21.annotate('',xy=(x1[int(n*f)],y1[int(n*f)]-.1),
                  xytext=(x1[int(n*f)],y1[int(n*f)]),
                  arrowprops=arrowdict)

    f=.95;n=len(x1)
    ax21.annotate('',xy=(x1[int(n*f)],y1[int(n*f)]+.4),
                  xytext=(x1[int(n*f)],y1[int(n*f)]),
                  arrowprops=arrowdict)
    
    ax21.annotate('',xy=(x1[int(n*f)],y1[int(n*f)]-.3),
                  xytext=(x1[int(n*f)],y1[int(n*f)]),
                  arrowprops=arrowdict)

    # plot arrows pointing towers llower stable manifold
    f=.34;n=len(x0)
    ax21.annotate('',xy=(x0[int(n*f)],y0[int(n*f)]-.02),
                  xytext=(x0[int(n*f)],-.4),
                  arrowprops=arrowdict)

    f=.66;
    ax21.annotate('',xy=(x0[int(n*f)],y0[int(n*f)]-.02),
                  xytext=(x0[int(n*f)],-.4),
                  arrowprops=arrowdict)


    # plot arrows pointing towers upper stable manifold
    f=.32;n=len(x2)
    ax21.annotate('',xy=(x2[int(n*f)],y2[int(n*f)]+.02),
                  xytext=(x2[int(n*f)],.4),
                  arrowprops=arrowdict)

    f=.68;
    ax21.annotate('',xy=(x2[int(n*f)],y2[int(n*f)]+.02),
                  xytext=(x2[int(n*f)],.4),
                  arrowprops=arrowdict)

    
    
    ax22_contour = ax22.contour(ze_rangem*ze_scale2,
                                Um*a.U_scale,D2,[0],colors='k')
    ax22.clear()
    
    p0=ax22_contour.collections[0].get_paths()[0]
    p1=ax22_contour.collections[0].get_paths()[1]
    p2=ax22_contour.collections[0].get_paths()[2]
    
    x0=p0.vertices[:,0];y0=p0.vertices[:,1]
    x1=p1.vertices[:,0];y1=p1.vertices[:,1]
    x2=p2.vertices[:,0];y2=p2.vertices[:,1]

    ax22.plot(x0,y0,color='k',lw=2)
    ax22.plot(x1,y1,color='tab:red',lw=2)
    ax22.plot(x2,y2,color='k',lw=2)


    
    ax23.plot(Z*a.Z_scale,ze_plot2*ze_scale2,color='black',lw=2)
    

    ax11.set_title(r'\textbf{A} Critical Manifold',x=0.2)
    ax11.set_xlabel(r'$Z$ (\si{\um})')
    ax11.set_ylabel(r'$U$ (\si{\um/s})')
    
    # make some labels scientific notation
    ax12.set_title(r'\textbf{B} Bifurcation Curve',x=0.2)
    ax12.set_xlabel(r'$\zeta$ (\si{kg/s})')
    ax12.set_ylabel(r'$U$ (\si{\um/s})')
    
    ax13.set_title(r'\textbf{C} Viscous Drag',x=.2)
    ax13.set_xlabel(r'$Z$ (\si{\um})')
    ax13.set_ylabel(r'$\zeta$ (\si{kg/s})')

    #https://matplotlib.org/api/_as_gen/...
    #matplotlib.axes.Axes.ticklabel_format.html
    ax13.ticklabel_format(axis='y',scilimits=(0,0),style='sci')

    ax21.set_title(r'\textbf{D} Critical Manifold',x=0.2)
    ax21.set_xlabel(r'$Z$ (\si{\um})')
    ax21.set_ylabel(r'$U$ (\si{\um/s})')
    
    
    ax22.set_title(r'\textbf{E} Bifurcation Curve',x=0.2)
    ax21.set_ylabel(r'$U$ (\si{\um/s})')
    
    ax22.set_xlabel(r'$\zeta$ (\si{kg/s})')
    
    ax23.set_title(r'\textbf{F} Viscous Drag',x=.2)
    ax23.set_xlabel(r'$Z$ (\si{\um})')
    ax23.set_ylabel(r'$\zeta$ (\si{kg/s})')
    
    leg = ax11.legend(prop={'size':8})
    leg.get_frame().set_linewidth(0.0)

    ax11.set_xlim(-5,5)
    ax12.set_xlim(0,20*ze_scale1)
    ax13.set_xlim(-5,5)

    ax21.set_xlim(-5,5)
    ax22.set_xlim(0,20*ze_scale2)
    ax23.set_xlim(-5,5)

    ax11.set_ylim(-.4,.4)
    ax12.set_ylim(-.4,.4)

    ax21.set_ylim(-.4,.4)
    ax22.set_ylim(-.4,.4)

    ax12.ticklabel_format(axis='x',scilimits=(0,0),style='sci')
    ax13.ticklabel_format(axis='y',scilimits=(0,0),style='sci')
    
    ax22.ticklabel_format(axis='x',scilimits=(0,0),style='sci')
    ax23.ticklabel_format(axis='y',scilimits=(0,0),style='sci')

    # needed to extract offset exponent
    # https://alexpearce.me/2014/04/exponent-label-in-matplotlib/    
    fig.savefig('junk.png') 
    
    offset = ax12.get_xaxis().get_offset_text()
    offset.set_visible(False)
    ax12.set_xlabel('{0} {1}'.format(ax12.get_xlabel(), offset.get_text()))

    offset = ax13.get_yaxis().get_offset_text()
    offset.set_visible(False)
    ax13.set_ylabel('{0} {1}'.format(ax13.get_ylabel(), offset.get_text()))

    offset = ax22.get_xaxis().get_offset_text()
    offset.set_visible(False)
    ax22.set_xlabel('{0} {1}'.format(ax22.get_xlabel(), offset.get_text()))

    offset = ax23.get_yaxis().get_offset_text()
    offset.set_visible(False)
    ax23.set_ylabel('{0} {1}'.format(ax23.get_ylabel(), offset.get_text()))

    plt.tight_layout()

    return fig




def add_fold_labels_ze(ax_ze1,ze1,dx=-.05,dy=0,
                       change_for_ze3=False,plot_idx=0):
    """
    add fold labels

    ax_ze1: axis object for ze slice.
    ze1: one-parameter bifurcation file
    dx,dy: offset for fold label

    plot_idx: just for convenience to make
    labeling multistability regions easier
    """
    
    # find all folds.
    B = ze1[:-1,0]-ze1[1:,0]
    
    #print(B[B==2])
    switch_idx1 = (B[:-2]==0)*(B[1:-1]==-1)*(B[2:]==0)
    switch_idx2 = (B[:-2]==0)*(B[1:-1]==1)*(B[2:]==0)

    switch1xs = ze1[1:-2][switch_idx1,3]
    switch1ys = ze1[1:-2][switch_idx1,6]

    switch2xs = ze1[1:-2][switch_idx2,3]
    switch2ys = ze1[1:-2][switch_idx2,6]

    print(switch1xs,switch1ys)
    print(switch2xs,switch2ys)

    ax_ze1.annotate("",xy=(switch1xs[0],switch1ys[0]),
                    xytext=(switch1xs[0]+dx,switch1ys[0]+dy),
                    arrowprops=dict(arrowstyle="->"))
    ax_ze1.annotate("",xy=(switch1xs[1],switch1ys[1]),
                    xytext=(switch1xs[1]+dx,switch1ys[1]+dy),
                    arrowprops=dict(arrowstyle="->"))
    ax_ze1.annotate("",xy=(switch2xs[0],switch2ys[0]),
                    xytext=(switch2xs[0]+dx,switch2ys[0]+dy),
                    arrowprops=dict(arrowstyle="->"))
    ax_ze1.annotate("",xy=(switch2xs[1],switch2ys[1]),
                    xytext=(switch2xs[1]+dx,switch2ys[1]+dy),
                    arrowprops=dict(arrowstyle="->"))

    
    if change_for_ze3:
        ax_ze1.scatter(switch1xs[1]+dx,switch1ys[1]+dy,marker='4',
                       color=purp,s=150,linewidth=3,zorder=-1)
        ax_ze1.scatter(switch1xs[0]+dx,switch1ys[0]+dy,marker='x',
                       color=lightpurp,s=150,linewidth=3,zorder=-1)
        ax_ze1.scatter(switch2xs[1]+dx,switch2ys[1]+dy,marker='*',
                       color=red,s=150,linewidth=3,zorder=-1)
        ax_ze1.scatter(switch2xs[0]+dx,switch2ys[0]+dy,marker='+',
                       color=yellow,s=150,linewidth=3,zorder=-1)

    

        
    else:
        ax_ze1.scatter(switch1xs[0]+dx,switch1ys[0]+dy,marker='4',
                       color=purp,s=150,linewidth=3,zorder=-1)
        ax_ze1.scatter(switch1xs[1]+dx,switch1ys[1]+dy,marker='x',
                       color=lightpurp,s=150,linewidth=3,zorder=-1)
        ax_ze1.scatter(switch2xs[0]+dx,switch2ys[0]+dy,marker='*',
                       color=red,s=150,linewidth=3,zorder=-1)
        ax_ze1.scatter(switch2xs[1]+dx,switch2ys[1]+dy,marker='+',
                       color=yellow,s=150,linewidth=3,zorder=-1)



    rect = []

    if plot_idx == 0:

        # switch1xs[0] # right, blue tristar
        # switch1xs[1] # left, blue X
        # switch2xs[0] # orange star
        # switch2xs[1] # +
        
        # show colors
        h = .2
        y_start = -.1
        
        # tristar -> right
        width1=.6-switch1xs[0];start1=(switch1xs[0],y_start) 
        
        # left -> X
        width2=.4-switch1xs[1];start2=(switch1xs[1],y_start) 
        
        # orange star -> tristar
        width3=switch1xs[0]-switch2xs[0];start3=(switch2xs[0],y_start) 
        
        # + -> orange star
        width4=switch2xs[0]-switch2xs[1];start4=(switch2xs[1],y_start)
        
        # X -> +
        width5=switch2xs[1]-switch1xs[1];start5=(switch1xs[1],y_start) 
        
        rect.append(patches.Rectangle(start1,width1,h,facecolor=color1pts,
                                      alpha=color1pts_al))
        rect.append(patches.Rectangle(start2,width2,h,facecolor=color1pts,
                                      alpha=color1pts_al))
        rect.append(patches.Rectangle(start3,width3,h,facecolor=color3pts,
                                      alpha=color3pts_al))
        rect.append(patches.Rectangle(start4,width4,h,facecolor=color1pts,
                                      alpha=color1pts_al))
        rect.append(patches.Rectangle(start5,width5,h,facecolor=color3pts,
                                      alpha=color3pts_al))

        col = PatchCollection(rect, zorder=-10,match_original=True)


    elif plot_idx == 1:

        # same as above
        # switch1xs[0] # right, blue tristar
        # switch1xs[1] # left, blue X
        # switch2xs[0] # orange star
        # switch2xs[1] # +
        
        # show colors
        h = .2
        y_start = -.1
        
        # tristar -> right
        width1=.6-switch1xs[0];start1=(switch1xs[0],y_start) 
        
        # left -> X
        width2=.4-switch1xs[1];start2=(switch1xs[1],y_start) 
        
        # X-> orange star
        width3=switch1xs[1]-switch2xs[0];start3=(switch2xs[0],y_start)
        
        # orange star -> +
        width4=switch2xs[0]-switch2xs[1];start4=(switch2xs[1],y_start)
        
        # + -> tristar
        width5=switch2xs[1]-switch1xs[0];start5=(switch1xs[0],y_start)
        
        rect.append(patches.Rectangle(start1,width1,h,facecolor=color1pts,
                                      alpha=color1pts_al))
        rect.append(patches.Rectangle(start2,width2,h,facecolor=color1pts,
                                      alpha=color1pts_al))
        rect.append(patches.Rectangle(start3,width3,h,facecolor=color3pts,
                                      alpha=color3pts_al))
        rect.append(patches.Rectangle(start4,width4,h,facecolor=color5pts,
                                      alpha=color5pts_al))
        rect.append(patches.Rectangle(start5,width5,h,facecolor=color3pts,
                                      alpha=color3pts_al))

        col = PatchCollection(rect, zorder=-10,match_original=True)


    elif plot_idx == 2:

        # switch1xs[0] # X
        # switch1xs[1] # tristar
        # switch2xs[0] # +
        # switch2xs[1] # orange star
        
        # show colors
        h = .2
        y_start = -.1
        
        # + -> right
        width1=0.6-switch2xs[0];start1=(switch2xs[0],y_start)
        
        # left -> orange star
        width2=0.4-switch2xs[1];start2=(switch2xs[1],y_start)
        
        # orange star -> X
        width3=switch2xs[1]-switch1xs[0];start3=(switch1xs[0],y_start)
        
        # X -> tristar
        width4=switch1xs[0]-switch1xs[1];start4=(switch1xs[1],y_start)
        
        # tristar -> +
        width5=switch1xs[1]-switch2xs[0];start5=(switch2xs[0],y_start) 
        
        rect.append(patches.Rectangle(start1,width1,h,facecolor=color1pts,
                                      alpha=color1pts_al))
        rect.append(patches.Rectangle(start2,width2,h,facecolor=color1pts,
                                      alpha=color1pts_al))
        rect.append(patches.Rectangle(start3,width3,h,facecolor=color3pts,
                                      alpha=color3pts_al))
        rect.append(patches.Rectangle(start4,width4,h,facecolor=color5pts,
                                      alpha=color5pts_al))
        rect.append(patches.Rectangle(start5,width5,h,facecolor=color3pts,
                                      alpha=color3pts_al))

        col = PatchCollection(rect, zorder=-10,match_original=True)


    ax_ze1.add_collection(col)


    return ax_ze1



def add_fold_labels_phi1(ax_phi1,phi1,dx=-.05,dy=0,change_for_ze3=False,
                         includelist1=[],includelist2=[],
                         markers1=[],markers2=[],
                         colors1=[],colors2=[],showall=False,plot_idx=0):
    
    """
    add fold labels
    """
    # find all folds.
    B = phi1[:-1,0]-phi1[1:,0]
    
    #print(B[B==2])
    switch_idx1 = (B[:-2]==0)*(B[1:-1]==-1)*(B[2:]==0)
    switch_idx2 = (B[:-2]==0)*(B[1:-1]==1)*(B[2:]==0)

    switch1xs = phi1[1:-2][np.roll(switch_idx1,1),3]
    switch1ys = phi1[1:-2][np.roll(switch_idx1,1),6]

    switch2xs = phi1[1:-2][np.roll(switch_idx2,1),3]
    switch2ys = phi1[1:-2][np.roll(switch_idx2,1),6]

    print('x1s',switch1xs,'y1s',switch1ys,'1idxs',
          np.arange(len(switch_idx1))[switch_idx1])
    print(switch2xs,switch2ys,np.arange(len(switch_idx2))[switch_idx2])

    # this is unreadable...
    j = 0
    for i in range(len(switch1xs)):
            
        # if showall is false, plot only those points specificed by
        # includelist. Includelist must be built manually.
        if not(showall): 
            if (i in includelist1):
                ax_phi1.scatter(switch1xs[i]+dx,switch1ys[i]+dy,
                                marker=markers1[j],color=colors1[j],
                                s=150,linewidth=3,zorder=-1)
                ax_phi1.annotate("",xy=(switch1xs[i],switch1ys[i]),
                                 xytext=(switch1xs[i]+dx,switch1ys[i]+dy),
                                 arrowprops=dict(arrowstyle="->"))

                j += 1
        else:
            ax_phi1.scatter(switch1xs[i]+dx,switch1ys[i])
            ax_phi1.annotate('on'+str(i),xy=(switch1xs[i],switch1ys[i]),
                             xytext=(switch1xs[i]+dx,switch1ys[i]+dy),
                             arrowprops=dict(arrowstyle="->"))
            


    j = 0
    for i in range(len(switch2xs)):

        if not(showall):
            if (i in includelist2):
                ax_phi1.scatter(switch2xs[i]+dx,switch2ys[i]+dy,
                                marker=markers2[j],color=colors2[j],
                                s=150,linewidth=3,zorder=-1)
                ax_phi1.annotate("",xy=(switch2xs[i],switch2ys[i]),
                                 xytext=(switch2xs[i]+dx,switch2ys[i]+dy),
                                 arrowprops=dict(arrowstyle="->"))
                j += 1

        else:
            ax_phi1.scatter(switch2xs[i]+dx,switch2ys[i]+dy)
            ax_phi1.annotate('off'+str(i),xy=(switch2xs[i],switch2ys[i]),
                             xytext=(switch2xs[i]+dx,switch2ys[i]+dy),
                             arrowprops=dict(arrowstyle="->"))

    rect = []

    if plot_idx == 0:

        p1 = [ switch1xs[i] for i in includelist1 ]
        p2 = [ switch2xs[i] for i in includelist2 ]


        # p1[0] +, p2[0] orange star, p2[1] X
        
        xmax = np.amax(phi1[:,3])
        
        h = .2
        y_start = -.1
        width1=0-p2[0];start1=(p2[0],y_start) # left -> orange star
        width2=p1[0]-xmax;start2=(xmax,y_start) # + -> right
        width3=p2[0]-p2[1];start3=(p2[1],y_start) # orange star -> X
        width4=p2[1]-p1[0];start4=(p1[0],y_start) # X -> +

        print('p1,p2',p1,p2,'xmax',xmax)
        
        rect.append(patches.Rectangle(start1,width1,h,facecolor=color3pts,
                                      alpha=color3pts_al))
        rect.append(patches.Rectangle(start2,width2,h,facecolor=color1pts,
                                      alpha=color1pts_al))
        rect.append(patches.Rectangle(start3,width3,h,facecolor=color1pts,
                                      alpha=color1pts_al))
        rect.append(patches.Rectangle(start4,width4,h,facecolor=color3pts,
                                      alpha=color3pts_al))
        
        col = PatchCollection(rect, zorder=-10,match_original=True)
    else:

        p1 = [ switch1xs[i] for i in includelist1 ]
        p2 = [ switch2xs[i] for i in includelist2 ]

        # p1[0] orange star +, p1[1] +
        
        xmax = np.amax(phi1[:,3])
        
        h = .2
        y_start = -.1
        width1=0-p1[0];start1=(p1[0],y_start) # left -> orange star
        width2=p1[0]-p1[1];start2=(p1[1],y_start) # orange star -> +
        width3=p1[1]-xmax;start3=(xmax,y_start) # + -> right
        

        rect.append(patches.Rectangle(start1,width1,h,facecolor=color5pts,
                                      alpha=color5pts_al))
        rect.append(patches.Rectangle(start2,width2,h,facecolor=color3pts,
                                      alpha=color3pts_al))
        rect.append(patches.Rectangle(start3,width3,h,facecolor=color1pts,
                                      alpha=color1pts_al))
        
        col = PatchCollection(rect, zorder=-10,match_original=True)
        
    ax_phi1.add_collection(col)


    return ax_phi1


def twopar_detailed():
    """
    the two parameter diagrams come from force.ode
    parameters:

    p p1=0.5,xi=.32
    p pi2=0.07195590236721254,pi3 = 1,pi4 = 4.7
    p pi5=0.1,pi6 = 10
    p2=1-p1

    the readme file in the same directory contains numerical info
    """
    
    
    a = lubrication()

    ze_scale = 1e-6*(6*np.pi*a.Rp*a.mu) # kg/s
    #z_scale = #6*np.pi*a.Rp*0.12
    print(ze_scale)

    #print(type_dict)

    fig = plt.figure(figsize=(8,7))
    gs = gridspec.GridSpec(3, 3)

    gs.update(hspace=.4, wspace=.5)
    ax_main = plt.subplot(gs[:-1,1:])
    ax_ze1 = plt.subplot(gs[0,0]) # first zeta slice
    ax_ze2 = plt.subplot(gs[1,0]) # second zeta slice
    ax_ze3 = plt.subplot(gs[2,0]) # third zeta slice

    ax_ph1 = plt.subplot(gs[-1,1])
    ax_ph2 = plt.subplot(gs[-1,2])

    twopar = np.loadtxt('bifurcations/allinfo_2par.dat')
    
    ########################## first zeta slice ze=7.83
    ze1 = np.loadtxt('bifurcations/ze=7.83.dat')

    #self.U_scale = self.F0/(6*pi*self.Rp*self.mu)
    
    #U_scale = oself.F0/(6*pi*self.Rp*self.mu)
    val_dict,type_dict = collect_disjoint_branches(ze1,sv_tol=.05,
                                                   redundant_threshold=0,
                                                   remove_isolated=True)
    keys = list(val_dict.keys())

    # dummy points to make legend easier.
    ax_ze1.plot([0,.5],[0,.5],color='k',label='Stable')
    ax_ze1.plot([0,.5],[0,.5],color='r',label='Unstable')

    ax_ze1 = add_fold_labels_ze(ax_ze1,ze1,dx=0,dy=-.03,plot_idx=0)
    
    for i in range(len(keys)):
        if np.sum(type_dict[keys[i]][1:,0]-1)==0:
            ls = '-'
            color = 'k'
            
        else:
            ls = '-'
            color = 'r'
        ax_ze1.plot(val_dict[keys[i]][:,0],val_dict[keys[i]][:,2],
                    color=color,lw=3,ls=ls)

    # label SN

    ########################## second zeta slice ze=5.86
    ze2 = np.loadtxt('bifurcations/ze=5.86.dat')
    val_dict,type_dict = collect_disjoint_branches(ze2,sv_tol=.05,
                                                   redundant_threshold=0,
                                                   remove_isolated=True)
    keys = list(val_dict.keys())

    for i in range(len(keys)):
        if np.sum(type_dict[keys[i]][1:,0]-1)==0:
            ls = '-'
            color = 'k'
        else:
            ls = '-'
            color = 'r'
        ax_ze2.plot(val_dict[keys[i]][:,0],val_dict[keys[i]][:,2],
                    color=color,label=keys[i],lw=3,ls=ls)

    ax_ze2 = add_fold_labels_ze(ax_ze2,ze2,dx=0,dy=-.05,plot_idx=1)

    ########################## third zeta slice ze=2.72
    ze3 = np.loadtxt('bifurcations/ze=2.72.dat.bk')
    val_dict,type_dict = collect_disjoint_branches(ze3,sv_tol=.05,
                                                   redundant_threshold=0,
                                                   remove_isolated=True)
    keys = list(val_dict.keys())

    for i in range(len(keys)):
        if np.sum(type_dict[keys[i]][1:,0]-1)==0:
            ls = '-'
            color = 'k'
        else:
            ls = '-'
            color = 'r'
        ax_ze3.plot(val_dict[keys[i]][:,0],val_dict[keys[i]][:,2],
                    color=color,label=keys[i],lw=3,ls=ls)

    ax_ze3 = add_fold_labels_ze(ax_ze3,ze3,dx=0,dy=-.025,
                                change_for_ze3=True,plot_idx=2)
        
    ########################## first phi slice phi=0.477
    ph1 = np.loadtxt('bifurcations/p1=.477.dat')
    ph1[:,3] *= ze_scale
    
    val_dict,type_dict = collect_disjoint_branches(ph1,sv_tol=1,
                                                   redundant_threshold=0,
                                                   remove_isolated=True)
    keys = list(val_dict.keys())

    for i in range(len(keys)):
        if np.sum(type_dict[keys[i]][1:,0]-1)==0:
            ls = '-'
            color = 'k'
        else:
            ls = '-'
            color = 'r'

        
        ax_ph1.plot(val_dict[keys[i]][:,0],val_dict[keys[i]][:,2],
                    color=color,label=keys[i],lw=3,ls=ls)    

    ax_ph1 = add_fold_labels_phi1(ax_ph1,ph1,dx=0,dy=-.05,
                                  includelist1=[2],markers1=['+'],
                                  colors1=[yellow],
                                  includelist2=[0,2],markers2=['*','x'],
                                  colors2=[red,lightpurp],
                                  plot_idx = 0)
        
    ########################## second phi slice phi=0.4912
    ph2 = np.loadtxt('bifurcations/p1=.4912.dat')
    ph2[:,3] *= ze_scale
    val_dict,type_dict = collect_disjoint_branches(ph2,sv_tol=1,
                                                   redundant_threshold=0,
                                                   remove_isolated=True)
    keys = list(val_dict.keys())

    for i in range(len(keys)):
        if np.sum(type_dict[keys[i]][1:,0]-1)==0:
            ls = '-'
            color = 'k'
        else:
            ls = '-'
            color = 'r'
        ax_ph2.plot(val_dict[keys[i]][:,0],val_dict[keys[i]][:,2],
                    color=color,lw=3,ls=ls)

    ax_ph2 = add_fold_labels_phi1(ax_ph2,ph2,dx=0,dy=-.038,
                                  includelist1=[0,2],markers1=['*','+'],
                                  colors1=[red,yellow],
                                  includelist2=[],markers2=[],
                                  colors2=[],showall=False,
                                  plot_idx = 1)

    ####################### two parameter diagram

    z1_value = 7.83*ze_scale
    z2_value = 5.86*ze_scale
    z3_value = 2.72*ze_scale

    phi_label = r'$\phi$' # r'$\phi$'
    ze_label = r'$\zeta$' # r'$\tilde\zeta$'
    U_label = r'$U$ (\si{\um/s})'
    
    # show slices for 1par bifurcations
    ax_main.plot([.477,.477],[0,10],ls='-',color='gray',
                 label=phi_label+' $=0.48$') # phi1 
    ax_main.plot([.4912,.4912],[0,10],ls='--',color='gray',
                 label=phi_label+' $=0.49$') # phi1 

    # split zi values into significant digit and exponents
    z1_value_sig,z1_value_exp = ('%.1E' % Decimal(str(z1_value))).split('E-')
    z2_value_sig,z2_value_exp = ('%.1E' % Decimal(str(z2_value))).split('E-')
    z3_value_sig,z3_value_exp = ('%.1E' % Decimal(str(z3_value))).split('E-')

    # z1_value_exp[1:] need to use this to ignore the leading 0.
    # e.g., 1e-06 becomes 1e-6
    z1_label = (r' $='+str(z1_value_sig)
                +r'\times 10^{-'+str(z1_value_exp[1:])+'}$')
    z2_label = (r' $='+str(z2_value_sig)
                +r'\times 10^{-'+str(z2_value_exp[1:])+'}$')
    z3_label = (r' $='+str(z3_value_sig)
                +r'\times 10^{-'+str(z3_value_exp[1:])+'}$')
    
    ax_main.plot([.4,.6],[z1_value,z1_value],ls='-.',
                 color='gray',label=ze_label+z1_label) # ze
    ax_main.plot([.4,.6],[z2_value,z2_value],ls=':',
                 color='gray',label=ze_label+z2_label) # ze
    ax_main.plot([.4,.6],[z3_value,z3_value],
                 dashes=[4, 1, 1, 1, 1, 1],
                 color='gray',label=ze_label+z3_label) # ze

    # number of points in each region
    ax_main.text(.5,2*ze_scale,'5',fontsize=20,horizontalalignment='center')
    ax_main.text(.55,2*ze_scale,'3',fontsize=20,horizontalalignment='center')
    ax_main.text(.45,2*ze_scale,'3',fontsize=20,horizontalalignment='center')
    ax_main.text(.55,5*ze_scale,'1',fontsize=20,horizontalalignment='center')
    ax_main.text(.45,5*ze_scale,'1',fontsize=20,horizontalalignment='center')
    ax_main.text(.5125,6.*ze_scale,'3',fontsize=20,
                 horizontalalignment='center')
    ax_main.text(1-.5125,6.*ze_scale,'3',fontsize=20,
                 horizontalalignment='center')
    ax_main.text(.5,9*ze_scale,'1',fontsize=20,horizontalalignment='center')
    
    val_dict,type_dict = collect_disjoint_branches(twopar,sv_tol=.05,
                                                   redundant_threshold=.00)
    
    keys = list(val_dict.keys())
    
    ax_main = fill_axis(ax_main,keys,val_dict,type_dict,ze_scale)

    ze_min = 0*ze_scale
    ze_max = 10*ze_scale
    
    phi_min = .4
    phi_max = .6

    u_min = -.07
    u_max = .07

    # match large figure to subfigures
    from matplotlib.patches import ConnectionPatch
    # ze1
    xy_ze1a = (phi_max,(u_max+u_min)/2)
    xy_ze1b = (phi_max,(u_max+u_min)/2)
    xy_main_ze1 = (phi_min,z1_value)

    con1 = ConnectionPatch(xyA=xy_ze1a,xyB=xy_main_ze1,
                           coordsA="data", coordsB="data",
                           axesA=ax_ze1,axesB=ax_main,color='gray')
    con2 = ConnectionPatch(xyA=xy_ze1b,xyB=xy_main_ze1,
                           coordsA="data", coordsB="data",
                           axesA=ax_ze1,axesB=ax_main,color='gray')

    ax_ze1.add_artist(con2)

    # ze 2
    xy_ze2a = (phi_max,(u_max+u_min)/2)
    xy_ze2b = (phi_max,(u_max+u_min)/2)
    xy_main_ze2 = (phi_min,z2_value)

    con1 = ConnectionPatch(xyA=xy_ze2a,xyB=xy_main_ze2, coordsA="data",
                           coordsB="data",
                           axesA=ax_ze2,axesB=ax_main,color='gray')
    con2 = ConnectionPatch(xyA=xy_ze2b,xyB=xy_main_ze2, coordsA="data",
                           coordsB="data",
                           axesA=ax_ze2,axesB=ax_main,color='gray')

    ax_ze2.add_artist(con2)


    # ze 3
    xy_ze3a = (phi_max,(u_max+u_min)/2)
    xy_ze3b = (phi_max,(u_max+u_min)/2)
    xy_main_ze3 = (phi_min,z3_value)

    con1 = ConnectionPatch(xyA=xy_ze3a,xyB=xy_main_ze3, coordsA="data",
                           coordsB="data",
                           axesA=ax_ze3,axesB=ax_main,color='gray')
    con2 = ConnectionPatch(xyA=xy_ze3b,xyB=xy_main_ze3, coordsA="data",
                           coordsB="data",
                           axesA=ax_ze3,axesB=ax_main,color='gray')

    ax_ze3.add_artist(con2)


    # phi1 1
    xy_p11a = ((ze_min+ze_max)/2,u_max)
    xy_p11b = ((ze_min+ze_max)/2,u_max)
    xy_main_p11 = (.477,ze_min)

    con1 = ConnectionPatch(xyA=xy_p11a,xyB=xy_main_p11, coordsA="data",
                           coordsB="data",
                           axesA=ax_ph1,axesB=ax_main,color='gray')
    con2 = ConnectionPatch(xyA=xy_p11b,xyB=xy_main_p11, coordsA="data",
                           coordsB="data",
                           axesA=ax_ph1,axesB=ax_main,color='gray')

    ax_ph1.add_artist(con2)


    # phi1 2
    xy_p12a = ((ze_min+ze_max)/2,u_max)
    xy_p12b = ((ze_min+ze_max)/2,u_max)
    xy_main_p12 = (.4912,ze_min)

    con1 = ConnectionPatch(xyA=xy_p12a,xyB=xy_main_p12, coordsA="data",
                           coordsB="data",
                           axesA=ax_ph2,axesB=ax_main,color='gray')
    con2 = ConnectionPatch(xyA=xy_p12b,xyB=xy_main_p12, coordsA="data",
                           coordsB="data",
                           axesA=ax_ph2,axesB=ax_main,color='gray')

    ax_ph2.add_artist(con2)

            
    # label saddle-nodes
    ax_main.text(.435,2*ze_scale,'SN')
    ax_main.text(.478,1*ze_scale,'SN')
    ax_main.text(.515,1*ze_scale,'SN')
    ax_main.text(.556,3*ze_scale,'SN')
            
    ax_main.set_xlabel(phi_label,fontsize=size)
    ax_main.set_ylabel(ze_label + ' (\si{kg/s})',fontsize=size)

    ax_ze1.legend()
    ax_main.legend(prop={'size':8})

    ax_ze1.set_xlim(phi_min,phi_max)
    ax_ze2.set_xlim(phi_min,phi_max)
    ax_ze3.set_xlim(phi_min,phi_max)
    ax_main.set_xlim(phi_min,phi_max)

    ax_ph1.set_xlim(0,10*ze_scale)
    ax_ph2.set_xlim(0,10*ze_scale)

    ax_ze1.set_ylim(u_min,u_max)
    ax_ze2.set_ylim(u_min,u_max)
    ax_ze3.set_ylim(u_min,u_max)
    ax_main.set_ylim(0,10*ze_scale)

    ax_ph1.set_ylim(u_min,u_max)
    ax_ph2.set_ylim(u_min,u_max)
    
    
    ax_ze1.set_title(r'\textbf{A} '+ze_label+z1_label
                     +' (\si{kg/s})',x=.2,fontsize=size)
    ax_ze2.set_title(r'\textbf{B} '+ze_label+z2_label
                     +' (\si{kg/s})',x=.2,fontsize=size)
    ax_ze3.set_title(r'\textbf{C} '+ze_label+z3_label
                     +' (\si{kg/s})',x=.2,fontsize=size)

    ax_main.set_title(r'\textbf{D}',x=.0,fontsize=size)
    
    ax_ph1.set_title(r'\textbf{E} '+phi_label+'$=0.48$',x=.2,fontsize=size)
    ax_ph2.set_title(r'\textbf{F} '+phi_label+'$=0.49$',x=.2,fontsize=size)

    #ax_ze1.set_xlabel('$\phi$')
    #ax_ze2.set_xlabel('$\phi_2$')
    ax_ze3.set_xlabel(phi_label,fontsize=size)

    ax_ph1.set_xlabel(ze_label + ' (\si{kg/s})',fontsize=size)
    ax_ph2.set_xlabel(ze_label + ' (\si{kg/s})',fontsize=size)

    ax_ze1.set_ylabel(U_label,fontsize=size)
    ax_ze2.set_ylabel(U_label,fontsize=size)
    ax_ze3.set_ylabel(U_label,fontsize=size)

    #ax_ph1.set_ylabel(r'$\tilde U$')
    #ax_ph2.set_ylabel(r'$\tilde U$')


    # force scientific notation
    ax_main.ticklabel_format(axis='y',scilimits=(0,0),style='sci')
    ax_ph1.ticklabel_format(axis='x',scilimits=(0,0),style='sci')
    ax_ph2.ticklabel_format(axis='x',scilimits=(0,0),style='sci')

    # needed to extract offset exponent
    # https://alexpearce.me/2014/04/exponent-label-in-matplotlib/
    fig.savefig('junk.png') 

    offset = ax_main.get_yaxis().get_offset_text()
    offset.set_visible(False)
    ax_main.set_ylabel('{0} {1}'.format(ax_main.get_ylabel(), offset.get_text()))
    
    offset = ax_ph1.get_xaxis().get_offset_text()
    offset.set_visible(False)
    ax_ph1.set_xlabel('{0} {1}'.format(ax_ph1.get_xlabel(), offset.get_text()))

    offset = ax_ph2.get_xaxis().get_offset_text()
    offset.set_visible(False)
    ax_ph2.set_xlabel('{0} {1}'.format(ax_ph2.get_xlabel(), offset.get_text()))

    plt.tight_layout()

    return fig




def pi4_vs_pi5():
    fig = plt.figure(figsize=(8,6))

    gs = gridspec.GridSpec(2,3)
    gs.update(hspace=.5, wspace=.5)

    ax11 = plt.subplot(gs[0,:2])
    ax12 = plt.subplot(gs[0,-1])
    ax21 = plt.subplot(gs[-1,0])

    ax22 = plt.subplot(gs[-1,1])
    ax23 = plt.subplot(gs[-1,2])
    
    pi3=1
    pi4=4.7
    pi5=.1
    
    pi6=10
    z=2

    ze_scale = 1e-6*(6*np.pi*0.96*1.2) # kg/s
    #ze_scale = 1 # kg/s
    
    p56=pi5/pi6
    ep4=exp(pi4)
    ep5=exp(pi5)
    
    ############################# main figure
    
    zvals = np.array([0,2,6,10])#np.linspace(2,10,4)
    pi4_array = np.linspace(0,5,50)
    pi5_array = np.zeros(len(pi4_array))

    ls = ['-','--','-.',':']

    for j in range(len(zvals)):

        if zvals[j] == 0:
            u = np.linspace(-.4,-.00,10000)
        else:
            u = np.linspace(-.04,-.00,1000)
            
        for i in range(len(pi4_array)):
            if zvals[j] == 0 and pi4_array[i]>=3.2:
                # there are some problems with root-finding when zeta = 0,
                # so cut off root finding for this pi4.
                pi5_array[i] = np.nan 
                
            else:
                try:
                    pi5_array[i] = brentq(zero_crossings,.01,.5,
                                          args=(u,zvals[j],pi3,
                                                pi4_array[i],pi6),rtol=1e-4)
                    
                except ValueError:
                    pi5_array[i] = np.nan


        # get scientific notation.
        print(('%.1E' % Decimal(str(zvals[j]*ze_scale))).split('E'))
        formatted = ('%.1E' % Decimal(str(zvals[j]*ze_scale)))
        ze_value_sig,ze_value_exp = formatted.split('E')
        if j == 0:
            label = r'$\zeta=0$'
        else:
            label = (r'$\zeta='+ze_value_sig
                     +r'\times 10^{-'+ze_value_exp[-1]+r'}$')
            
        ax11.plot(pi4_array,pi5_array,label=label,ls=ls[j],color='k')

    # label subplot locations
    dx=.1;dy=-.01
    ax11.scatter(3, 0.15, color='k',s=10)
    ax11.scatter(1.5, 0.3, color='k',s=10)
    
    ax11.text(3+dx, 0.15+dy, r'\textbf{B}')
    ax11.text(1.5+dx, 0.3+dy, r'\textbf{C}')
    
    ax11.scatter(3, 0.3, color='k',s=10)
    ax11.scatter(1.5, 0.15, color='k',s=10)
    ax11.scatter(4.7, 0.1, color='k',s=80,marker='*')
    
    ax11.text(3+dx, 0.3+dy, r'\textbf{D}')
    ax11.text(1.5+dx, 0.15+dy, r'\textbf{E}')

    ax11.set_xlabel(r'$\pi_4$',fontsize=12)
    ax11.set_ylabel(r'$\pi_5$',fontsize=12)
    
    ax11.set_title(r'\textbf{A} Cusp Bifurcations')
    ax11.legend(loc='upper right',ncol=2)
    
    ax11.set_xlim(0,5)
    ax11.set_ylim(0,.55)
    
    
    ze_min = 0
    ze_max = 6*ze_scale
        
    ################### pi4=3, pi5=0.15
    dat = np.loadtxt('bifurcations/pi4=3_pi5=0.15.dat')
    
    x = dat[:,7][dat[:,2]==0]
    y = dat[:,3][dat[:,2]==0]

    skipn = 6
    ax12.scatter(x[::skipn],y[::skipn]*ze_scale,s=1,color='k')

    # get max zeta and show go/no-go
    ax12.set_xlim(0.4,0.6)
    ax12.set_ylim(ze_min,ze_max)
    ax12.set_ylabel(r'$\zeta$  (\si{kg/s})',fontsize=12)
    #ax12.set_xlabel('$\phi$')

    ax12.set_title(r'\textbf{B} $\pi_4=3$, $\pi_5=0.15$')

    
    ################### pi4=1.5, pi5=0.3
    dat = np.loadtxt('bifurcations/pi4=1.5_pi5=0.3.dat')

    x = dat[:,7][dat[:,2]==0]
    y = dat[:,3][dat[:,2]==0]

    skipn = 2
    ax21.scatter(x[::skipn],y[::skipn]*ze_scale,s=1,color='k')

    ax21.set_xlim(0.4,0.6)
    ax21.set_ylim(-0.1*ze_scale,ze_max)
    ax21.set_xlabel('$\phi$',fontsize=12)

    ax21.set_title(r'\textbf{C} $\pi_4=1.5$, $\pi_5=0.3$')
    ax21.set_ylabel(r'$\zeta$ (\si{kg/s})',fontsize=12)

    
    ################### pi4=3, pi5=0.3
    dat = np.loadtxt('bifurcations/pi4=3_pi5=0.3.dat')

    x = dat[:,7][dat[:,2]==0]
    y = dat[:,3][dat[:,2]==0]

    skipn = 3
    ax22.scatter(x[::skipn],y[::skipn]*ze_scale,s=1,color='k')

    ax22.set_xlim(0.4,0.6)
    ax22.set_ylim(ze_min,ze_max)
    ax22.set_xlabel('$\phi$',fontsize=12)

    ax22.set_ylabel(r'$\zeta$  (\si{kg/s})',fontsize=12)

    ax22.set_title(r'\textbf{D} $\pi_4=3$, $\pi_5=0.3$')

    
    ################### pi4=1.5, pi5=0.15
    dat = np.loadtxt('bifurcations/pi4=1.5_pi5=0.15_v2.dat')

    x = dat[:,3][dat[:,2]==1]
    y = dat[:,4][dat[:,2]==1]

    skipn = 3
    ax23.scatter(x[::skipn],y[::skipn]*ze_scale,s=1,color='k')

    ax23.set_ylabel(r'$\zeta$  (\si{kg/s})',fontsize=12)

    ax12.set_xlim(0.4,0.6)
    ax21.set_xlim(0.4,0.6)
    ax22.set_xlim(0.4,0.6)
    ax23.set_xlim(0.4,0.6)

    ax23.set_ylim(ze_min,ze_max)
    
    ax23.set_xlabel('$\phi$',fontsize=12)
    ax23.set_title(r'\textbf{E} $\pi_4=1.5$, $\pi_5=0.15$')

    ax11.tick_params(axis='both',which='major', labelsize=12)
    ax12.tick_params(axis='both',which='major', labelsize=12)
    ax21.tick_params(axis='both',which='major', labelsize=12)
    ax22.tick_params(axis='both',which='major', labelsize=12)
    ax23.tick_params(axis='both',which='major', labelsize=12)

    ax12.ticklabel_format(axis='y',scilimits=(0,0),style='sci')
    ax21.ticklabel_format(axis='y',scilimits=(0,0),style='sci')
    ax22.ticklabel_format(axis='y',scilimits=(0,0),style='sci')
    ax23.ticklabel_format(axis='y',scilimits=(0,0),style='sci')
    
    # needed to extract offset exponent
    # https://alexpearce.me/2014/04/exponent-label-in-matplotlib/
    fig.savefig('junk.png') 
    
    offset = ax12.get_yaxis().get_offset_text()
    offset.set_visible(False)
    ax12.set_ylabel('{0} {1}'.format(ax12.get_ylabel(), offset.get_text()))

    offset = ax21.get_yaxis().get_offset_text()
    offset.set_visible(False)
    ax21.set_ylabel('{0} {1}'.format(ax21.get_ylabel(), offset.get_text()))

    offset = ax22.get_yaxis().get_offset_text()
    offset.set_visible(False)
    ax22.set_ylabel('{0} {1}'.format(ax22.get_ylabel(), offset.get_text()))

    offset = ax23.get_yaxis().get_offset_text()
    offset.set_visible(False)
    ax23.set_ylabel('{0} {1}'.format(ax23.get_ylabel(), offset.get_text()))

    
    return fig



def minimal_2par():
    """
    minimal 2 parameter figure for use in final figure in discussion.
    """

    a = lubrication()
    ze_scale = 1e-6*(6*np.pi*a.Rp*a.mu) # kg/s

    fig = plt.figure(figsize=(3,3))
    ax = fig.add_subplot(111)
    
    twopar = np.loadtxt('bifurcations/allinfo_2par.dat')

    val_dict,type_dict = collect_disjoint_branches(twopar,sv_tol=.05,
                                                   redundant_threshold=.00)
    keys = list(val_dict.keys())

    ax = fill_axis(ax,keys,val_dict,type_dict,ze_scale,tol=1e-3,beautify=False)
    ax.set_yscale('log')
    
    ax.set_ylabel(r'$\zeta$ (\si{kg/s})')
    ax.set_xlabel(r'$\phi$')

    ax.set_ylim(3e-6,2e-2)
    ax.set_xlim(.4,.6)

    return fig


def cylinder():
    
    """
    TODO: add small cartoon axes to subplots.
    mmake dome head for spine figure
    
    """
    
    T1 = .1
    
    gs = gridspec.GridSpec(nrows=2,ncols=3,wspace=-.1,hspace=.5)
    fig = plt.figure(figsize=(5,4))
    ax11 = fig.add_subplot(gs[:,:2],projection='3d')
    ax12 = fig.add_subplot(gs[0,2])
    ax22 = fig.add_subplot(gs[1,2])
    
    
    a = lubrication(phi1=.57,Rp=.96,Rc=1.22,
                    pi3=1,pi4=4.7,pi5=0.1,pi6=10,
                    mu=1.2,T=T1,constriction='piecewise',U0=0.2,
                    dt=0.02,eps=1,
                    F0=50,method='euler')
    a.Z0 = -5/a.Rp
    
    z = np.linspace(-7,7,100)  # dimensional
    r = a.pi1(z)
    th = np.linspace(0,2*np.pi,100)
    
    radius_al = 0.25
    
    # draw arrow going into spine
    
    ar1 = Arrow3D([0,0],[0,0],[-5,-1],
                  mutation_scale=10, 
                  lw=2, arrowstyle="-|>", color="k")
    
    ax11.add_artist(ar1)
    
    # draw spine
    Z,TH = np.meshgrid(z,th)
    #Z,TH = np.mgrid[-7:7:.1, 0:2*np.pi:.1]
    X = np.zeros_like(Z)
    Y = np.zeros_like(Z)
    print(np.shape(Z))
    for i in range(len(Z[:,0])):
        X[i,:] = a.pi1(Z[i,:])*np.cos(TH[i,:])
        Y[i,:] = a.pi1(Z[i,:])*np.sin(TH[i,:])
    
    ax11.plot_surface(X,Y,Z,alpha=.25)
    
    shifts = np.array([3,-3,0])
    names = ['x','y','z']
    size = 2
    
    
    for i in range(3):
        coords = np.zeros((3,2))
        
        coords[:,0] += shifts
        coords[:,1] += shifts
        
        coords[i][1] += size
        arx = Arrow3D(*list(coords),
                      mutation_scale=5, 
                      lw=2, arrowstyle="-|>", color="k")
    
        ax11.text(*list(coords[:,1]),names[i],horizontalalignment='center')
    
        ax11.add_artist(arx)
        
    

    # draw sphere for cap
    b = a.base_radius
    r = np.sqrt(b**2+7**2)
    th2 = np.linspace(0,np.arctan(b/7),100)
    phi = np.linspace(0,2*np.pi,100)
    
    TH2,PHI = np.meshgrid(th2,phi)
    X = r*np.sin(TH2)*np.cos(PHI)
    Y = r*np.sin(TH2)*np.sin(PHI)
    Z = r*np.cos(TH2)
    ax11.plot_surface(X,Y,Z,color='tab:blue',alpha=.5)

 
    # highlight radius change
    ax11.plot(a.base_radius*np.cos(th),
              a.base_radius*np.sin(th),
              7,color='tab:blue',alpha=radius_al)
    ax11.plot(a.base_radius*np.cos(th),
              a.base_radius*np.sin(th),
              -5,color='tab:blue',alpha=radius_al)
    ax11.plot(a.base_radius*np.cos(th),
              a.base_radius*np.sin(th),
              5,color='tab:blue',alpha=radius_al)
    
    ax11.plot(a.Rc*np.cos(th),
              a.Rc*np.sin(th),
              -a.inner_width/2,color='tab:blue',alpha=radius_al)
    ax11.plot(a.Rc*np.cos(th),
              a.Rc*np.sin(th),
              a.inner_width/2,color='tab:blue',alpha=radius_al)
    
    # draw sphere vesicle
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    X = np.cos(u)*np.sin(v)
    Y = np.sin(u)*np.sin(v)
    Z = np.cos(v)
    ax11.plot_surface(X,Y,Z,color='gray',alpha=.5)
    
    # label spine head and base
    ax11.text(5,-5,6,'Spine Head')
    
    
    ax11.text(5,-5,-6,'Spine Base')
    
    
    # draw plane 1
    plane1_color = 'green'
    X,Y = np.meshgrid(np.linspace(-5,5,2),np.linspace(-5,5,2))
    ax11.plot_surface(X,Y,np.zeros_like(X),color=plane1_color,alpha=0.0)
    
    # draw plane 2 (coronal)
    #ax11.plot_surface(X,Y,Z,color='k')
    plane1_color = 'green'
    Y,Z = np.meshgrid(np.linspace(-7,7,10),np.linspace(-7,7,10))
    ax11.plot_surface(np.zeros_like(Y),Y,Z,color=plane1_color,alpha=0.0)
    
    
    # B
    # Draw vesicle and spine wall
    ax12.plot(a.Rc*np.cos(th),a.Rc*np.sin(th),color='tab:blue',alpha=.5)
    ax12.fill(a.Rp*np.cos(th),a.Rp*np.sin(th),color='gray',alpha=.5)
    
    ax12.text(-.6,.4,'Vesicle')
    ax12.text(.9,.9,'Spine Wall')
    
    # annotate Rp length
    ax12.annotate('', xy=(0, -a.Rc), xytext=(0, 0), xycoords='data',
               arrowprops=dict(arrowstyle='-|>',color='k', lw=1))
    ax12.text(-.35,-.5,r'$R_c$')

    # annotate Rc length channel
    ax12.annotate(r'', xy=(a.Rp, 0.0), xytext=(0, 0), xycoords='data', 
                  arrowprops=dict(arrowstyle='-|>', color='k',lw=1))
    ax12.text(.4,.2,r'$R_p$')
    
    # annotate Rc length channel
    ax12.annotate(r'$h(z)$', xy=(1.08,-.1), xytext=(1.08, -1.5),
                  xycoords='data',fontsize=10, ha='center', va='center',
                  arrowprops=dict(arrowstyle=('-[, widthB='+str(a.Rc-a.Rp+.2)
                                              +', lengthB=.5'), lw=.5))

    # draw axes
    shift = .75
    x_center, y_center = (-1.5,-1.5)

    ax12.annotate(r'x', xy=(x_center,y_center),
                  xytext=(x_center+shift,y_center),
                  va='center',
                  arrowprops=dict(mutation_scale=5, 
                                  arrowstyle='<|-', 
                                  color='k',lw=2),annotation_clip=False)
    
    ax12.annotate(r'y', xy=(x_center,y_center),
                  xytext=(x_center,y_center+shift),
                  ha='center', 
                  arrowprops=dict(mutation_scale=5,
                                  arrowstyle='<|-',
                                  color='k',lw=2),annotation_clip=False)
    
    
    
    # C
    # draw molecular motors
    
    pad = .15
    ax22.fill(a.Rp*np.cos(th),a.Rp*np.sin(th),color='gray',alpha=.5)
    
    ax22.plot([a.Rc+pad,a.Rc+pad],[-2,2],color='tab:blue',alpha=.5)
    ax22.plot([-a.Rc-pad,-a.Rc-pad],[-2,2],color='tab:blue',alpha=.5)
    
    
    ax22.text(0,0,'Vesicle',ha='center',va='center')
    ax22.text(-a.Rc-3*pad,0,'Spine Wall',rotation=90)
    
    ax22.annotate(r'Center of Mass $Z$',xy=(a.Rp,0),xytext=(a.Rc+4*pad,0),
                  ha='center',va='center',rotation=-90,
                  arrowprops=dict(arrowstyle='-|>',
                                  color='k',lw=1),annotation_clip=False)
    
    # draw axes
    shift = .75
    x_center, y_center = (-1.5,-1.5)

    ax22.annotate(r'x,y', xy=(x_center,y_center),
                  xytext=(x_center+shift,y_center),
                  va='center',
                  arrowprops=dict(mutation_scale=5, 
                                  arrowstyle='<|-', 
                                  color='k',lw=2),
                  annotation_clip=False)
    
    ax22.annotate(r'z', xy=(x_center,y_center),
                  xytext=(x_center,y_center+shift),
                  ha='center', 
                  arrowprops=dict(mutation_scale=5,
                                  arrowstyle='<|-',
                                  color='k',lw=2),
                  annotation_clip=False)
    
    
    ax11.set_title(r'\textbf{A} Idealized Spine Geometry',x=.45,y=1.09)
    ax12.set_title(r'\textbf{B} Transverse Cross-section',x=.5)
    ax22.set_title(r'\textbf{C} Molecular Motors',x=.27)
    
    
    # set equal aspect ratios
    ax12.set_aspect(1)
    ax22.set_aspect(1)
    
    ax11.set_axis_off()
    ax12.set_axis_off()
    ax22.set_axis_off()
    
    ax22.set_xticks([])
    ax22.set_yticks([])
    
    lo = -4.4
    hi = 4.4
    
    dx = -.5
    
    ax11.set_xlim(lo-dx,hi+dx)
    ax11.set_ylim(lo-dx,hi+dx)
    ax11.set_zlim(lo,hi)
    
    ax22.set_xlim(-1.4,1.4)
    ax22.set_ylim(-1.4-pad,1.4+pad)
    
    ax11.view_init(20,45)
    
    return fig

def fv():
    gs = gridspec.GridSpec(nrows=2,ncols=3,wspace=.5,hspace=.7)
    fig = plt.figure(figsize=(8,4))
    ax11 = fig.add_subplot(gs[0,0])
    ax12 = fig.add_subplot(gs[0,1])
    ax13 = fig.add_subplot(gs[0,2])

    ax21 = fig.add_subplot(gs[1,0])
    ax22 = fig.add_subplot(gs[1,1])
    ax23 = fig.add_subplot(gs[1,2])
    
    
    # parset 1 
    # A,D
    kwargs = {'phi1':.51,'Rp':.96,'Rc':1.22,
              'pi3':1,'pi4':4.7,'pi5':0.1,'pi6':10,
              'mu':1.2,'T':.01,'constriction':'piecewise','U0':0.2,
              'dt':0.02,'eps':1,'ze':6,'F0':50,'method':'euler'}
    
    a = lubrication(**kwargs)
    a.Z0 = -5/a.Rp
    
    Us = np.linspace(-.1,.1,5000)
    
    title = r'\textbf{A} $\phi='+str(a.phi1)+'$, $\zeta='+str(a.ze)+'$'
    ax11.set_title(title,x=0.4)
    ax11.plot(Us,a.phi1*a.FAm(Us),color='tab:blue',
              label='$\phi F_{-A}(U)$',ls='--',dashes=[2,1])
    ax11.plot(Us,(1-a.phi1)*a.FA(Us),color='tab:red',
              label='$(1-\phi)F_A(U)$',ls='--',dashes=[2,1])
    ax11.plot(Us,a.ze*Us,color='gray',ls='-',lw=1,
              label=r'$\zeta U$')
    ax11.plot(Us,a.F(Us),color='tab:purple',
              label=r'$\phi F_{-A}(U) + (1-\phi)F_{A}(U)$')
    
    # total force
    total_f = a.F(Us)-a.ze*Us

    # get zeros
    c_idxs = np.where(np.diff(np.sign(total_f)))[0][1:-1]
    ax21.plot([Us[0],Us[-1]],[0,0],color='gray',lw=1)
    ax21.plot(Us,a.F(Us)-a.ze*Us,color='tab:green')
    
    # label crossings    
    ax21.scatter(Us[c_idxs],np.zeros(len(Us[c_idxs])),color='k',
                 zorder=6,s=6,label='Force-Balance')
    ax11.scatter(Us[c_idxs],a.F(Us)[c_idxs],color='k',zorder=6,s=6)
    
    # draw stability arrows
    ax21.arrow(-.03,0,.02,0,
               head_width=0.02,head_length=0.005,zorder=9,
               fc='k')
    ax21.arrow(.045,0,-.005,0,
               head_width=0.02,head_length=0.005,zorder=6,
               fc='k')
    
    # inset
    axins = inset_axes(ax21, width="60%", height="30%", loc=3)
    axins.plot([Us[0],Us[-1]],[0,0],color='gray',lw=1)
    axins.plot(Us,(a.F(Us)-a.ze*Us),color='tab:green')
    axins.scatter(Us[c_idxs],np.zeros(len(Us[c_idxs])),color='k',zorder=6,s=6)    
    mark_inset(ax21, axins, loc1=2, loc2=4, fc="none", ec='0.5',alpha=0.5)
    # draw stability arrows
    axins.arrow(-.005,0,.001,0,
               head_width=0.005,head_length=0.002,zorder=9,
               fc='k')
    axins.arrow(.006,0,-.001,0,
               head_width=0.005,head_length=0.002,zorder=6,
               fc='k')
    axins.arrow(.015,0,.001,0,
               head_width=0.005,head_length=0.002,zorder=6,
               fc='k')
    axins.arrow(.033,0,-.001,0,
               head_width=0.005,head_length=0.002,zorder=6,
               fc='k')
    axins.set_xticks([])
    axins.set_yticks([])
    axins.set_xlim(-.01,.035)
    axins.set_ylim(-.02,.02)
    
    
    # parset 2
    # B,E
    kwargs['phi1'] = 0.57
    kwargs['ze'] = 8
    a.__init__(**kwargs)
    
    Us = np.linspace(-.1,.1,1000)
    title = r'\textbf{B} $\phi='+str(a.phi1)+'$, $\zeta='+str(a.ze)+'$'
    ax12.set_title(title,x=0.4)
    ax12.plot(Us,a.phi1*a.FAm(Us),color='tab:blue',
              label='$\phi F_{-A}(U)$',ls='--',dashes=[2,1])
    ax12.plot(Us,(1-a.phi1)*a.FA(Us),color='tab:red',
              label='(1-\phi)$F_A(U)$',ls='--',dashes=[2,1])
    ax12.plot(Us,a.ze*Us,color='gray',lw=1,
              label=r'$\zeta U$')
    ax12.plot(Us,a.F(Us),color='tab:purple',
              label=r'$\phi F_{-A}(U) + (1-\phi)F_{A}(U)$')
    
    # total force
    total_f = a.F(Us)-a.ze*Us

    # get zeros
    c_idxs = np.where(np.diff(np.sign(total_f)))[0][1:-1]
    ax22.plot([Us[0],Us[-1]],[0,0],color='gray',lw=1)
    ax22.plot(Us,a.F(Us)-a.ze*Us,color='tab:green')
    
    # label crossings    
    ax22.scatter(Us[c_idxs],np.zeros(len(Us[c_idxs])),color='k',zorder=6,s=6)
    ax12.scatter(Us[c_idxs],a.F(Us)[c_idxs],color='k',zorder=6,s=6)
    
    
    # draw stability arrows
    ax22.arrow(-.025,0,.03,0,
               head_width=0.04,head_length=0.01,zorder=9,
               fc='k')
    ax22.arrow(.085,0,-.03,0,
               head_width=0.04,head_length=0.01,zorder=6,
               fc='k')
    
    # parset 3
    # C,F
    kwargs['phi1'] = .48
    kwargs['ze'] = 3
    a.__init__(**kwargs)
    
    Us = np.linspace(-.1,.1,1000)
    
    title = r'\textbf{C} $\phi='+str(a.phi1)+'$, $\zeta='+str(a.ze)+'$'
    ax13.set_title(title,x=0.4)
    ax13.plot(Us,a.phi1*a.FAm(Us),color='tab:blue',
              label='$\phi F_{-A}(U)$',ls='--',dashes=[2,1])
    ax13.plot(Us,(1-a.phi1)*a.FA(Us),color='tab:red',
              label='(1-\phi)$F_A(U)$',ls='--',dashes=[2,1])
    ax13.plot(Us,a.ze*Us,color='gray',lw=1,
              label=r'$\zeta U$')
    ax13.plot(Us,a.F(Us),color='tab:purple',
              label=r'$\phi F_{-A}(U) + (1-\phi)F_{A}(U)$')
    
    
    
    ax23.plot([Us[0],Us[-1]],[0,0],color='gray',lw=1)
    
    # total force
    total_f = a.F(Us)-a.ze*Us

    # get zeros
    c_idxs = np.where(np.diff(np.sign(total_f)))[0][1:-1]
    ax23.plot(Us,total_f,color='tab:green')
    
    # label crossings    
    ax23.scatter(Us[c_idxs],np.zeros(len(Us[c_idxs])),color='k',zorder=6,s=6)
    ax13.scatter(Us[c_idxs],a.F(Us)[c_idxs],color='k',zorder=6,s=6)
    
    # draw stability arrows
    ax23.arrow(-.095,0,.01,0,
               head_width=0.02,head_length=0.01,zorder=9,
               fc='k')
    ax23.arrow(-.02,0,-.02,0,
               head_width=0.02,head_length=0.01,zorder=6,
               fc='k')
    ax23.arrow(.02,0,.01,0,
               head_width=0.02,head_length=0.01,zorder=9,
               fc='k')
    ax23.arrow(.093,0,-.02,0,
               head_width=0.02,head_length=0.01,zorder=6,
               fc='k')
    
    # end par3
    
    # tweaks and labels
    ax11.legend(bbox_to_anchor=(1.9, 1.5), loc='upper center',ncol=4)
    ax21.legend(loc='upper right',fontsize=8)
    
    
    ax11.set_xlabel('U (Nondimensional)')
    ax11.set_ylabel('F (Nondimensional)')
    
    ax12.set_xlabel('U (Nondimensional)')
    ax12.set_ylabel('F (Nondimensional)')
    
    ax13.set_xlabel('U (Nondimensional)')
    ax13.set_ylabel('F (Nondimensional)')
    
    ax21.set_xlabel('U (Nondimensional)')
    ax21.set_ylabel('F (Nondimensional)')
    
    ax22.set_xlabel('U (Nondimensional)')
    ax22.set_ylabel('F (Nondimensional)')
    
    ax23.set_xlabel('U (Nondimensional)')
    ax23.set_ylabel('F (Nondimensional)')
    
    #ax11.set_xlim(Us[0],Us[-1])
    ax11.set_xlim(-.05,.05)
    ax12.set_xlim(Us[0],Us[-1])
    ax13.set_xlim(Us[0],Us[-1])
    
    ax21.set_xlim(-.05,.05)
    ax22.set_xlim(Us[0],Us[-1])
    ax23.set_xlim(Us[0],Us[-1])
    
    ax11.set_ylim(-.55,.55)
    
    ax21.set_ylim(-.2,.2)
    ax22.set_ylim(-.4,.4)
    ax23.set_ylim(-.2,.2)
    
    
    ax21.set_title(r'\textbf{D} Total Force',x=0.3)
    ax22.set_title(r'\textbf{E} Total Force',x=0.3)
    ax23.set_title(r'\textbf{F} Total Force',x=0.3)
    #ax11.plot(Us,a.F(Us))
    
    
    return fig

def generate_figure(function, args, filenames, dpi=100):
    # workaround for python bug where forked processes use the same random 
    # filename.
    #tempfile._name_sequence = None;

    fig = function(*args)

    if type(filenames) == list:
        for name in filenames:
            if name.split('.')[-1] == 'ps':
                fig.savefig(name, orientation='landscape',dpi=dpi)
            else:
                fig.savefig(name,dpi=dpi,bbox_inches='tight')
    else:
        if name.split('.')[-1] == 'ps':
            fig.savefig(filenames,orientation='landscape',dpi=dpi)
        else:
            fig.savefig(filenames,dpi=dpi)


def main():
    
    # listed in order of Figures in paper
    figures = [
        (cylinder,[],['cylinder.pdf']),
        (fv,[],['fv.pdf']),
        (constriction,[],["constriction.pdf"]),
        (critical_manifold_with_ze,[],["critical_manifold.pdf"]),
        (twopar_detailed,[],["bifurcations.pdf","bifurcations.png"]), 
        (pi4_vs_pi5,[],["pi4_vs_pi5.pdf"]),
        (minimal_2par,[],["minimal_2par.pdf"]),
    ]
    
    for fig in figures:
        generate_figure(*fig)


if __name__ == "__main__":
    main()
