# -*- coding: utf-8 -*-
"""
Created in 2024 by Carlos Forsythe

v1.0
-3 fitting types: uniform grid, shift rows, shift columns
-3 depths of search
-solution is centered, and minimum diameter given

 
"""

import tkinter as tk
import numpy as np
from scipy import optimize
import math
import time
import matplotlib.pyplot as plt
plt.ioff()
from matplotlib.patches import Rectangle, Circle
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

#%% defalut values
width_default=33
height_default=26
xspacing_default=.1
yspacing_default=.1
waferdiameter_default=100
edgeexclusionwidth_default=1

fittype_default=3 # 0-indexed from top
searchdepth_default=0

#%% other options
developer_plot=False #adds extra plotting elements

no_local_search_shift=True #no need for local search when output is integer

maxiter_grid0=1000
maxiter_grid1=2000
maxiter_grid2=4000

maxiter_shift0=1000
maxiter_shift1=10000
maxiter_shift2=60000

note='Note:\n'+\
    'This program uses global optimzers to find potential\n'+\
    'best solutions. This process is random and might\n'+\
    'return a different number of dies when repeated. This\n'+\
    'is generally only an issue when the best solution is\n'+\
    'a very tight fit. Increasing search depth can help ensure\n'+\
    'an optimal fit but Quick is generally recommended.'
#%% create window and frames
root = tk.Tk() # Create the root window
root.wm_attributes('-topmost', True) #put window on top
root.title('Dies Per Wafer v1.0') # Set window title
root.geometry('1200x900') # Set window size
#root.configure(bg='white')

frame1=tk.Frame(root)
#frame1.configure(bg='white')
frame1.grid(column = 1, row = 1,padx=50)
frame2=tk.Frame(root)
#frame1.configure(bg='white')
frame2.grid(column = 2, row = 1)

#%% make plot axes
fig, ax = plt.subplots(figsize=(8,8))
canvas= FigureCanvasTkAgg(fig,master=frame2)
canvas.get_tk_widget().pack()
toolbar = NavigationToolbar2Tk (canvas,frame2,pack_toolbar=False)
toolbar.update()
toolbar.pack(anchor='w',fill=tk.X)
#%% input numbers   
label_inputdirections = tk.Label(frame1,text = 'Input Dimensions (use uniform units)')
label_inputdirections.grid(column = 1, row = 1,pady=10, columnspan=2)
 
label_width = tk.Label(frame1,text = 'die width')
label_width.grid(column = 1, row = 2,sticky='E', ipadx=10)
width_s=tk.StringVar(value=str(width_default))
entry_width=tk.Entry(frame1,textvariable=width_s)
entry_width.grid(column = 2, row = 2)

label_height = tk.Label(frame1,text = 'die height')
label_height.grid(column = 1, row = 3,sticky='E', ipadx=10)
height_s=tk.StringVar(value=str(height_default))
entry_height=tk.Entry(frame1,textvariable=height_s)
entry_height.grid(column = 2, row = 3)

label_xspacing = tk.Label(frame1,text = 'x spacing')
label_xspacing.grid(column = 1, row = 4,sticky='E', ipadx=10)
xspacing_s=tk.StringVar(value=str(xspacing_default))
entry_xspacing=tk.Entry(frame1,textvariable=xspacing_s)
entry_xspacing.grid(column = 2, row = 4)

label_yspacing = tk.Label(frame1,text = 'y spacing')
label_yspacing.grid(column = 1, row = 5,sticky='E', ipadx=10)
yspacing_s=tk.StringVar(value=str(yspacing_default))
entry_yspacing=tk.Entry(frame1,textvariable=yspacing_s)
entry_yspacing.grid(column = 2, row = 5)

label_waferdiameter = tk.Label(frame1,text = 'wafer diameter')
label_waferdiameter.grid(column = 1, row = 6,sticky='E', ipadx=10)
waferdiameter_s=tk.StringVar(value=str(waferdiameter_default))
entry_waferdiameter=tk.Entry(frame1,textvariable=waferdiameter_s)
entry_waferdiameter.grid(column = 2, row = 6)

label_edgeexclusionwidth = tk.Label(frame1,text = 'edge exclusion width')
label_edgeexclusionwidth.grid(column = 1, row = 7,sticky='E', ipadx=10)
edgeexclusionwidth_s=tk.StringVar(value=str(edgeexclusionwidth_default))
entry_edgeexclusionwidth=tk.Entry(frame1,textvariable=edgeexclusionwidth_s)
entry_edgeexclusionwidth.grid(column = 2, row = 7)

#%% calculate optimized fit
spacer1 = tk.Label(frame1, text='')
spacer1.grid(column=1, row=8)

label_fittype = tk.Label(frame1,text = 'Fit Type:')
label_fittype.grid(column = 1, row = 9, columnspan=2,sticky='W', ipadx=60)
fittype_var = tk.IntVar(value=fittype_default)
FT0 = tk.Radiobutton(frame1, text='Uniform Grid', variable=fittype_var, value=0)
FT0.grid(column = 1, row = 10, columnspan=2,sticky='W', ipadx=60)
FT1 = tk.Radiobutton(frame1, text='Shift Rows', variable=fittype_var, value=1)
FT1.grid(column = 1, row = 11, columnspan=2,sticky='W', ipadx=60)
FT2 = tk.Radiobutton(frame1, text='Shift Columns', variable=fittype_var, value=2)
FT2.grid(column = 1, row = 12, columnspan=2,sticky='W', ipadx=60)
FT3 = tk.Radiobutton(frame1, text='all types (show best)', variable=fittype_var, value=3)
FT3.grid(column = 1, row = 13, columnspan=2,sticky='W', ipadx=60)

spacer2 = tk.Label(frame1, text='')
spacer2.grid(column=1, row=14)

label_fittype = tk.Label(frame1,text = 'Search Depth:')
label_fittype.grid(column = 1, row = 15, columnspan=2,sticky='W', ipadx=60)
searchdepth_var = tk.IntVar(value=searchdepth_default)
MI0 = tk.Radiobutton(frame1, text='Quick', variable=searchdepth_var, value=0)
MI0.grid(column = 1, row = 16, columnspan=2,sticky='W', ipadx=60)
MI1 = tk.Radiobutton(frame1, text='Thorough', variable=searchdepth_var, value=1)
MI1.grid(column = 1, row = 17, columnspan=2,sticky='W', ipadx=60)
MI2 = tk.Radiobutton(frame1, text='Exhaustive', variable=searchdepth_var, value=2)
MI2.grid(column = 1, row = 18, columnspan=2,sticky='W', ipadx=60)

def execute():
    button_execute.config(relief='sunken')
    extraoutput.set('searching...')
    root.update()
    print('\nnew fit...')
    width=float(width_s.get())
    height=float(height_s.get())
    xspacing=float(xspacing_s.get())
    yspacing=float(yspacing_s.get())
    waferdiameter=float(waferdiameter_s.get())
    edgeexclusionwidth=float(edgeexclusionwidth_s.get())
    fittype=fittype_var.get() #rewritten later if using all fits
    searchdepth=searchdepth_var.get()
    if searchdepth==0:
        maxiter_grid=maxiter_grid0
        maxiter_shift=maxiter_shift0
    elif searchdepth==1:
        maxiter_grid=maxiter_grid1
        maxiter_shift=maxiter_shift1
    elif searchdepth==2:
        maxiter_grid=maxiter_grid2
        maxiter_shift=maxiter_shift2
    ewr=waferdiameter/2-edgeexclusionwidth
    Nx=math.ceil(ewr/width) #spacing ignored
    Ny=math.ceil(ewr/height)
    Xoff,Yoff=np.meshgrid(range(-Nx,Nx+1),range(-Ny,Ny+1),indexing='ij')
    
    def calcXYcen(offsets,ft):
        if ft==0:
            return (offsets[0]+Xoff*(width+xspacing),
                    offsets[1]+Yoff*(height+yspacing))
        elif ft==1:
            return (np.tile(offsets[1:], (2*Nx+1,1))+Xoff*(width+xspacing),
                    offsets[0]+Yoff*(height+yspacing))
        elif ft==2:
            return (offsets[0]+Xoff*(width+xspacing),
                    np.tile(offsets[1:], (2*Ny+1,1)).transpose()+Yoff*(height+yspacing))
        
    def constructV(offsets,ft):
        Xcen,Ycen=calcXYcen(offsets,ft)
        return np.stack(((Xcen-width/2)**2+(Ycen-height/2)**2,
                         (Xcen-width/2)**2+(Ycen+height/2)**2,
                         (Xcen+width/2)**2+(Ycen-height/2)**2,
                         (Xcen+width/2)**2+(Ycen+height/2)**2))**0.5
    
    def withpartialscore(offsets):
        V=constructV(offsets,0)
        valid=np.all(V<=ewr,axis=0)
        PS=1-np.sum(np.clip(V-ewr,0,np.inf), axis=0)/4/(width**2+height**2)**0.5
        partial_score=max(PS[~valid])
        return -(valid.sum()+partial_score)
    
    def trueoffsets(reducedoffset,ft):
        if ft==0:
            return reducedoffset
        elif ft==1:
            return np.concatenate((reducedoffset[:1],(reducedoffset[1:]>0.5)*(width+xspacing)/2))
        elif ft==2:
            return np.concatenate((reducedoffset[:1],(reducedoffset[1:]>0.5)*(height+yspacing)/2))
        
    def fullonly_rowshift(reducedoffsets):
        V=constructV(trueoffsets(reducedoffsets,1),1)
        valid=np.all(V<=ewr,axis=0)
        return -valid.sum()
    
    def fullonly_colshift(reducedoffsets):
        V=constructV(trueoffsets(reducedoffsets,2),2)
        valid=np.all(V<=ewr,axis=0)
        return -valid.sum()
    
    start = time.time()
    if fittype==0:
        fit=optimize.dual_annealing(withpartialscore, [(0,(width+xspacing)/2), (0, (height+yspacing)/2)],
                                    maxiter=maxiter_grid) #, minimizer_kwargs={'tol': 1e-10}
        print('Uniform Grid fit output:')
        print(fit)
        
    elif fittype==1:
        fit= optimize.dual_annealing(fullonly_rowshift,
                                             [(0, (height+yspacing)/2)]+[(0,1) for _ in range(-Ny,Ny+1)],
                                             maxiter=maxiter_shift,no_local_search=no_local_search_shift)
        print('Row Shift fit output:')
        print(fit)
    elif fittype==2:
        fit= optimize.dual_annealing(fullonly_colshift,
                                             [(0, (width+xspacing)/2)]+[(0,1) for _ in range(-Nx,Nx+1)],
                                             maxiter=maxiter_shift,no_local_search=no_local_search_shift)
        print('Column Shift fit output:')
        print(fit)
    elif fittype==3:
        fit0=optimize.dual_annealing(withpartialscore, [(0,(width+xspacing)/2), (0, (height+yspacing)/2)],
                                     maxiter=maxiter_grid)
        print('Uniform Grid fit output:')
        print(fit0)
        fit1= optimize.dual_annealing(fullonly_rowshift,
                                             [(0, (height+yspacing)/2)]+[(0,1) for _ in range(-Ny,Ny+1)],
                                             maxiter=maxiter_shift,no_local_search=no_local_search_shift)
        print('\nRow Shift fit output:')
        print(fit1)
        fit2= optimize.dual_annealing(fullonly_colshift,
                                             [(0, (width+xspacing)/2)]+[(0,1) for _ in range(-Nx,Nx+1)],
                                             maxiter=maxiter_shift,no_local_search=no_local_search_shift)
        print('\nColumn Shift fit output:')
        print(fit2)
        N0=math.floor(-fit0.fun)
        N1=math.floor(-fit1.fun)
        N2=math.floor(-fit2.fun)
        if N0 == max((N0,N1,N2)):
            fittype=0
            fit=fit0
        elif N1 == max((N0,N1,N2)):
            fittype=1
            fit=fit1
        elif N2 == max((N0,N1,N2)):
            fittype=2
            fit=fit2
    end = time.time()
    Nfit=math.floor(-fit.fun)
    
    fit_offsets=trueoffsets(fit.x,fittype)
    fit_V=constructV(fit_offsets,fittype)
    fit_valid=np.tile(np.all(fit_V<=ewr,axis=0),[4,1,1])
    offsets_test=fit_offsets
    def Rmax(offsets,ft,validmask):
        if fittype==0:
            V=constructV(offsets,ft)
        else:
            offsets_test[0]=offsets #just 1 offset
            V=constructV(offsets_test,ft)
        return V[validmask].max()
    
    print('\nCentering output:')
    if fittype==0:
        fit2 = optimize.minimize(Rmax,fit_offsets,bounds=[(0,(width+xspacing)/2), (0, (height+yspacing)/2)],
                                 args=(fittype,fit_valid))
        centered_offsets=fit2.x
    elif fittype==1:
        fit2 = optimize.minimize_scalar(Rmax,bounds=(0, (height+yspacing)/2),args=(fittype,fit_valid))
        centered_offsets=fit_offsets
        centered_offsets[0]=fit2.x
    elif fittype==2:
        fit2 = optimize.minimize_scalar(Rmax,bounds=(0, (width+xspacing)/2),args=(fittype,fit_valid))
        centered_offsets=fit_offsets
        centered_offsets[0]=fit2.x
    print(fit2)
    
    min_diameter=2*(fit2.fun+edgeexclusionwidth)
    final_offsets=centered_offsets
    print('\nFinal Offsets:')
    print(final_offsets)
    print('Solution diameter:')
    print(min_diameter)
    area_utilization = Nfit*width*height*4*waferdiameter**-2/np.pi
    extraoutput.set('area utilization = {:.2%}\nminimum wafer diameter = {:.1f}\nsearch time: {:.2f} sec'\
                    .format(area_utilization,math.ceil(10*min_diameter)/10,end-start))

    Xcen,Ycen=calcXYcen(final_offsets,fittype)
    V=constructV(final_offsets,fittype)
    valid=np.all(V<=ewr,axis=0)
    partial=np.logical_and(np.any(V<=ewr,axis=0),np.logical_not(valid))
    ax.clear()
    if fittype==0:
        ax.set_title('Uniform Grid - '+str(Nfit)+' dies')
    elif fittype==1:
        ax.set_title('Shift Rows - '+str(Nfit)+' dies')
    elif fittype==2:
        ax.set_title('Shift Columns - '+str(Nfit)+' dies')
    for idx,_ in np.ndenumerate(Xcen):
        if partial[idx]:
            ax.add_patch(Rectangle((Xcen[idx]-width/2, Ycen[idx]-height/2),width, height,
                 edgecolor = 'xkcd:tomato', lw=1,fill=0))
    for idx,_ in np.ndenumerate(Xcen):
        if valid[idx]:
            ax.add_patch(Rectangle((Xcen[idx]-width/2, Ycen[idx]-height/2),width, height,
                 edgecolor = 'xkcd:blue', lw=1,facecolor='xkcd:light blue'))
    ax.add_patch(Circle((0,0),waferdiameter/2,edgecolor = 'k',fill=0))
    ax.add_patch(Circle((0,0),waferdiameter/2-edgeexclusionwidth,edgecolor = 'k',ls='--',fill=0))
    
    #optional plot elements
    if developer_plot:
        ax.add_patch(Circle((0,0),fit2.fun,edgecolor = 'xkcd:blue',ls='--',fill=0))
        ax.add_patch(Rectangle((0, 0),(width+xspacing)/2, (height+yspacing)/2,
             edgecolor = 'k',ls=':', lw=1,fill=0))
        if fittype==0:
            plt.scatter(final_offsets[0],final_offsets[1],marker='+',color='k')
        elif fittype==1:
            x=final_offsets[1:]
            y=final_offsets[0]+(height+yspacing)*Yoff[0,:]
            plt.scatter(x,y,marker='+',color='k')
        elif fittype==2:
            y=final_offsets[1:]
            x=final_offsets[0]+(width+xspacing)*Xoff[:,0]
            plt.scatter(x,y,marker='+',color='k')

    ax.set_aspect('equal')
    ax.autoscale()
    canvas.draw()
    button_execute.config(relief='raised')

button_execute=tk.Button(frame1,text='Find Best Fit',command=execute)
button_execute.grid(column = 1, row = 19,pady=30,ipadx=20,ipady=10, columnspan=2)

extraoutput= tk.StringVar()
label_extraoutput = tk.Label(frame1,textvariable = extraoutput)
label_extraoutput.grid(column = 1, row = 20, columnspan=2)

label_note = tk.Label(frame1,text = note)
label_note.grid(column = 1, row = 21,columnspan=2,rowspan=7,pady=20)

#%% excute window
root.mainloop()