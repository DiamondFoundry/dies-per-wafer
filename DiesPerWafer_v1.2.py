# -*- coding: utf-8 -*-
"""
Created in 2024 by Carlos Forsythe

v1.0
-3 fitting types: uniform grid, shift rows, shift columns
-3 depths of search
-solution is centered, and minimum diameter given

v1.1
-added shift&rot fit
-made fits checkboxes

v1.2
-added symmetry constraint
-allowed for blank height and yspacing
    
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
width_default=20
height_default=''
xspacing_default=1.0
yspacing_default=''
waferdiameter_default=200
edgeexclusionwidth_default=1

searchdepth_default=0

#select default fit checks with True/False
check_grid_default=True
check_ShiftRows_default=True
check_ShiftCols_default=True
check_ShiftRot_default=False

check_symmetric_default=False
#%% other options
allow_rotation=True

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
#root.wm_attributes('-topmost', True) #put window on top, annoying?
root.title('Dies Per Wafer v1.2') # Set window title
root.geometry('1200x900') # Set window size

frame1=tk.Frame(root)
frame1.grid(column = 1, row = 1,padx=50)
frame2=tk.Frame(root)
frame2.grid(column = 2, row = 1)

#%% make plot axes
#fig, ax = plt.subplots(figsize=(8,8))
fig=plt.figure(figsize=(8,8))
ax=fig.add_axes((.05,.05,.9,.85))
canvas= FigureCanvasTkAgg(fig,master=frame2)
canvas.get_tk_widget().pack()
toolbar = NavigationToolbar2Tk (canvas,frame2,pack_toolbar=False)
toolbar.update()
toolbar.pack(anchor='w',fill=tk.X)
#%% input numbers   
r=1
label_inputdirections = tk.Label(frame1,text = 'Input Dimensions (use uniform units)')
label_inputdirections.grid(column = 1, row = r,pady=10, columnspan=2)

r+=1
label_width = tk.Label(frame1,text = 'die width')
label_width.grid(column = 1, row = r,sticky='E', ipadx=10)
width_s=tk.StringVar(value=str(width_default))
entry_width=tk.Entry(frame1,textvariable=width_s)
entry_width.grid(column = 2, row = r)

r+=1
label_height = tk.Label(frame1,text = 'die height*')
label_height.grid(column = 1, row = r,sticky='E', ipadx=10)
height_s=tk.StringVar(value=str(height_default))
entry_height=tk.Entry(frame1,textvariable=height_s)
entry_height.grid(column = 2, row = r)

r+=1
label_xspacing = tk.Label(frame1,text = 'x spacing')
label_xspacing.grid(column = 1, row = r,sticky='E', ipadx=10)
xspacing_s=tk.StringVar(value=str(xspacing_default))
entry_xspacing=tk.Entry(frame1,textvariable=xspacing_s)
entry_xspacing.grid(column = 2, row = r)

r+=1
label_yspacing = tk.Label(frame1,text = 'y spacing*')
label_yspacing.grid(column = 1, row = r,sticky='E', ipadx=10)
yspacing_s=tk.StringVar(value=str(yspacing_default))
entry_yspacing=tk.Entry(frame1,textvariable=yspacing_s)
entry_yspacing.grid(column = 2, row = r)

r+=1
label_waferdiameter = tk.Label(frame1,text = 'wafer diameter')
label_waferdiameter.grid(column = 1, row = r,sticky='E', ipadx=10)
waferdiameter_s=tk.StringVar(value=str(waferdiameter_default))
entry_waferdiameter=tk.Entry(frame1,textvariable=waferdiameter_s)
entry_waferdiameter.grid(column = 2, row = r)

r+=1
label_edgeexclusionwidth = tk.Label(frame1,text = 'edge exclusion width')
label_edgeexclusionwidth.grid(column = 1, row = r,sticky='E', ipadx=10)
edgeexclusionwidth_s=tk.StringVar(value=str(edgeexclusionwidth_default))
entry_edgeexclusionwidth=tk.Entry(frame1,textvariable=edgeexclusionwidth_s)
entry_edgeexclusionwidth.grid(column = 2, row = r)

r+=1
starnote = tk.Label(frame1,text = '* leave empty to inherit above value')
starnote.grid(column = 1, row = r,columnspan=2, ipadx=10)

#%% make fit parameter selections
r+=1
spacer1 = tk.Label(frame1, text='')
spacer1.grid(column=1, row=r)
r+=1
label_fittype = tk.Label(frame1,text = 'Fit Types:')
label_fittype.grid(column = 1, row = r, columnspan=2,sticky='W', ipadx=60)
r+=1
boolvar_grid=tk.BooleanVar(value=check_grid_default)
check_grid = tk.Checkbutton(frame1, text='Uniform Grid',variable=boolvar_grid, onvalue=True, offvalue=False)
check_grid.grid(column = 1, row = r, columnspan=2,sticky='W', ipadx=80)
r+=1
boolvar_ShiftRows=tk.BooleanVar(value=check_ShiftRows_default)
check_ShiftRows = tk.Checkbutton(frame1, text='Shift Rows',variable=boolvar_ShiftRows, onvalue=True, offvalue=False)
check_ShiftRows.grid(column = 1, row = r, columnspan=2,sticky='W', ipadx=80)
r+=1
boolvar_ShiftCols=tk.BooleanVar(value=check_ShiftCols_default)
check_ShiftCols = tk.Checkbutton(frame1, text='Shift Columns',variable=boolvar_ShiftCols, onvalue=True, offvalue=False)
check_ShiftCols.grid(column = 1, row = r, columnspan=2,sticky='W', ipadx=80)
r+=1
boolvar_ShiftRot=tk.BooleanVar(value=check_ShiftRot_default)
check_ShiftRot = tk.Checkbutton(frame1, text='Shift & Rotate Rows',variable=boolvar_ShiftRot, onvalue=True, offvalue=False)
check_ShiftRot.grid(column = 1, row = r, columnspan=2,sticky='W', ipadx=80)
r+=1
spacer2 = tk.Label(frame1, text='')
spacer2.grid(column=1, row=r)
r+=1
boolvar_symmetry=tk.BooleanVar(value=check_symmetric_default)
check_symmetry = tk.Checkbutton(frame1, text='X & Y Mirror Symmetry',variable=boolvar_symmetry, onvalue=True, offvalue=False)
check_symmetry.grid(column = 1, row = r, columnspan=2,sticky='W', ipadx=60)
r+=1
spacer2 = tk.Label(frame1, text='')
spacer2.grid(column=1, row=r)
r+=1
label_fittype = tk.Label(frame1,text = 'Search Depth:')
label_fittype.grid(column = 1, row = r, columnspan=2,sticky='W', ipadx=60)
searchdepth_var = tk.IntVar(value=searchdepth_default)
MI0 = tk.Radiobutton(frame1, text='Quick', variable=searchdepth_var, value=0)
r+=1
MI0.grid(column = 1, row = r, columnspan=2,sticky='W', ipadx=80)
MI1 = tk.Radiobutton(frame1, text='Thorough', variable=searchdepth_var, value=1)
r+=1
MI1.grid(column = 1, row = r, columnspan=2,sticky='W', ipadx=80)
MI2 = tk.Radiobutton(frame1, text='Exhaustive', variable=searchdepth_var, value=2)
r+=1
MI2.grid(column = 1, row = r, columnspan=2,sticky='W', ipadx=80)

#%% Gather user inputs (on button press)
def execute():
    button_execute.config(relief='sunken')
    extraoutput.set('searching...')
    root.update()
    print('\nnew fit...')
    width=float(width_s.get())
    try:
        height=float(height_s.get())
    except ValueError:
        height=width
    xspacing=float(xspacing_s.get())
    try:
        yspacing=float(yspacing_s.get())
    except ValueError:
        yspacing=xspacing
    waferdiameter=float(waferdiameter_s.get())
    edgeexclusionwidth=float(edgeexclusionwidth_s.get())
    ft_grid=boolvar_grid.get()
    ft_ShiftRows=boolvar_ShiftRows.get()
    ft_ShiftCols=boolvar_ShiftCols.get()
    ft_ShiftRot=boolvar_ShiftRot.get()
    symmetric=boolvar_symmetry.get()
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
#%% find best solution
    Nx=math.ceil(ewr/width) #spacing ignored
    Ny=math.ceil(ewr/height)
    Nmax=max((Nx,Ny))
    Xoff,Yoff=np.meshgrid(range(-Nx,Nx+1),range(-Ny,Ny+1),indexing='ij') #grids of values [-Nx,...,Nx],[-Ny,...,Ny]
    Xoff2,Yoff2=np.meshgrid(range(-Nmax,Nmax+1),range(-Nmax,Nmax+1),indexing='ij') #grid for shift&rotate

    def CalculatePositions(offsets,ft): #return Xcorner,Ycorner, Xextent, Yextent 
        if ft==0:
            if symmetric:
                return ((offsets[0]>0.5)*(width+xspacing)/2+Xoff*(width+xspacing)-width/2,
                        (offsets[1]>0.5)*(height+yspacing)/2+Yoff*(height+yspacing)-height/2,
                        width*np.ones((2*Nx+1,2*Ny+1)),
                        height*np.ones((2*Nx+1,2*Ny+1)))
            else:
                return (offsets[0]+Xoff*(width+xspacing)-width/2,
                        offsets[1]+Yoff*(height+yspacing)-height/2,
                        width*np.ones((2*Nx+1,2*Ny+1)),
                        height*np.ones((2*Nx+1,2*Ny+1)))
        elif ft==1:
            if symmetric:
                if offsets[:1]<0.5: #center row unique
                    mirror_offset=np.concatenate((np.flip(offsets[2:]),offsets[1:]))
                else: #center row copied
                    mirror_offset=np.concatenate((np.flip(offsets[1:-1]),offsets[1:]))
                xoffsets=(mirror_offset>0.5)*(width+xspacing)/2
                yoffset=(offsets[:1]>0.5)*(height+yspacing)/2
            else:
                xoffsets=(offsets[1:]>0.5)*(width+xspacing)/2
                yoffset=offsets[0]
            return (np.tile(xoffsets, (2*Nx+1,1))+Xoff*(width+xspacing)-width/2,
                    yoffset+Yoff*(height+yspacing)-height/2,
                    width*np.ones((2*Nx+1,2*Ny+1)),
                    height*np.ones((2*Nx+1,2*Ny+1)))
        elif ft==2:
            if symmetric:
                if offsets[:1]<0.5: #center col unique
                    mirror_offset=np.concatenate((np.flip(offsets[2:]),offsets[1:]))
                else: #center col copied
                    mirror_offset=np.concatenate((np.flip(offsets[1:-1]),offsets[1:]))
                yoffsets=(mirror_offset>0.5)*(height+yspacing)/2
                xoffset=(offsets[:1]>0.5)*(width+xspacing)/2
            else:
                yoffsets=(offsets[1:]>0.5)*(height+yspacing)/2
                xoffset=offsets[0]                
            return (xoffset+Xoff*(width+xspacing)-width/2,
                    np.tile(yoffsets, (2*Ny+1,1)).transpose()+Yoff*(height+yspacing)-height/2,
                    width*np.ones((2*Nx+1,2*Ny+1)),
                    height*np.ones((2*Nx+1,2*Ny+1)))   
        elif ft==3:
            if symmetric:
                if offsets[:1]<0.5: #center row unique
                    xshifted_bool=np.concatenate((np.flip(offsets[2:Nmax+2]),offsets[1:Nmax+2]))>0.5
                    rotated_bool=np.concatenate((np.flip(offsets[Nmax+3:]),offsets[Nmax+2:]))>0.5
                else: #center row copied
                    xshifted_bool=np.concatenate((np.flip(offsets[1:Nmax+1]),offsets[1:Nmax+2]))>0.5
                    rotated_bool=np.concatenate((np.flip(offsets[Nmax+2:-1]),offsets[Nmax+2:]))>0.5
                yoffset0=(offsets[0]>0.5)*((height*~rotated_bool[Nmax]+width*rotated_bool[Nmax])+yspacing)/2
            else:
                xshifted_bool=offsets[1:2*Nmax+2]>0.5
                rotated_bool=offsets[2*Nmax+2:]>0.5
                yoffset0=offsets[0]
            xextent=~rotated_bool*width+rotated_bool*height #1D
            yextent=rotated_bool*width+~rotated_bool*height #1D
            xoffset0=xshifted_bool*xspacing/2-~xshifted_bool*xextent/2 #1D
            Rsum=np.zeros(2*Nmax+1) #1D
            Rsum[Nmax+1:]=rotated_bool[Nmax:-1].cumsum()
            Rsum[:Nmax]=-np.flip(np.flip(rotated_bool[:Nmax]).cumsum())
            notRsum=np.arange(-Nmax,Nmax+1)-Rsum
            yoffsets=yoffset0+yspacing*np.arange(-Nmax,Nmax+1)-\
                (height*~rotated_bool[Nmax]+width*rotated_bool[Nmax])/2+\
                Rsum*width+notRsum*height
            R=np.tile(rotated_bool, (2*Nmax+1,1))
            return(np.tile(xoffset0, (2*Nmax+1,1))+Xoff2*(~R*width+R*height+xspacing),
                   np.tile(yoffsets, (2*Nmax+1,1)),
                   np.tile(xextent,(2*Nmax+1,1)),
                   np.tile(yextent,(2*Nmax+1,1)))   
    def constructV2(offsets,ft): #return distance^2 of each corner from center of wafer
        Xcorner,Ycorner, Xextent, Yextent = CalculatePositions(offsets,ft)
        return np.stack((Xcorner**2+Ycorner**2,
                         (Xcorner+Xextent)**2+Ycorner**2,
                         Xcorner**2+(Ycorner+Yextent)**2,
                         (Xcorner+Xextent)**2+(Ycorner+Yextent)**2))
    def gridwithpartialscore(offsets,ft):
        V=constructV2(offsets,ft)
        valid=np.all(V<=ewr**2,axis=0)
        PS=1-np.sum(np.clip(V**0.5-ewr,0,np.inf), axis=0)/4/(width**2+height**2)**0.5
        partial_score=max(PS[~valid])
        return -(valid.sum()+partial_score)
    def countfulldies(offsets,ft):
        V=constructV2(offsets,ft)
        valid=np.all(V<=ewr**2,axis=0)
        return -valid.sum()
    
    if all(~np.array([ft_grid,ft_ShiftRows,ft_ShiftCols,ft_ShiftRot])):
        ax.clear()
        canvas.draw()
        extraoutput.set('No Fit Types Selected')
        button_execute.config(relief='raised')
        return
    
    Nfits=np.zeros(4)
    start = time.time()
    start2 = start
    if ft_grid:
        if symmetric:
            #there are literally only 4 solutions in this case, so this is kinda dumb
            fit0=optimize.dual_annealing(countfulldies, [(0,1), (0,1)],
                                        args=(0,),maxiter=maxiter_grid)
        else:
            fit0=optimize.dual_annealing(gridwithpartialscore, [(0,(width+xspacing)/2), (0, (height+yspacing)/2)],
                                        args=(0,),maxiter=maxiter_grid)
        print('Uniform Grid fit output:')
        print(fit0)
        print('Uniform Grid time: {:.2f} sec'.format(time.time()-start2))
        Nfits[0]=math.floor(-fit0.fun)
        start2=time.time()
    if ft_ShiftRows:
        if symmetric:
            fit1= optimize.dual_annealing(countfulldies,
                                                 [(0, 1)]+[(0,1) for _ in range(0,Ny+1)],
                                                 args=(1,),maxiter=maxiter_shift,no_local_search=no_local_search_shift)
        else:
            fit1= optimize.dual_annealing(countfulldies,
                                                 [(0, (height+yspacing)/2)]+[(0,1) for _ in range(-Ny,Ny+1)],
                                                 args=(1,),maxiter=maxiter_shift,no_local_search=no_local_search_shift)    
        print('Shift Rows fit output:')
        print(fit1)
        print('Shift Rows time: {:.2f} sec'.format(time.time()-start2))
        Nfits[1]=math.floor(-fit1.fun)
        start2=time.time()
    if ft_ShiftCols:
        if symmetric:
            fit2= optimize.dual_annealing(countfulldies,
                                                 [(0, 1)]+[(0,1) for _ in range(0,Nx+1)],
                                                 args=(2,),maxiter=maxiter_shift,no_local_search=no_local_search_shift)
        else:
            fit2= optimize.dual_annealing(countfulldies,
                                                 [(0, (width+xspacing)/2)]+[(0,1) for _ in range(-Nx,Nx+1)],
                                                 args=(2,),maxiter=maxiter_shift,no_local_search=no_local_search_shift)
        print('Shift Columns fit output:')
        print(fit2)
        print('Shift Columns time: {:.2f} sec'.format(time.time()-start2))
        Nfits[2]=math.floor(-fit2.fun)
        start2=time.time()
    if ft_ShiftRot:
        if symmetric:
            fit3= optimize.dual_annealing(countfulldies,
                                                 [(0, 1)]+[(0,1) for _ in range(2*Nmax+2)],
                                                 args=(3,),maxiter=maxiter_shift,no_local_search=no_local_search_shift)
        else:
            fit3= optimize.dual_annealing(countfulldies,
                                                 [(0, (max((height,width))+xspacing)/2)]+[(0,1) for _ in range(4*Nmax+2)],
                                                 args=(3,),maxiter=maxiter_shift,no_local_search=no_local_search_shift)
        print('Shifted & Rotated fit output:')
        print(fit3)
        print('Shifted & Rotated time: {:.2f} sec'.format(time.time()-start2))
        Nfits[3]=math.floor(-fit3.fun)
        start2=time.time()
        
    fittype=Nfits.argmax() #"indices corresponding to the first occurrence are returned"
    if fittype==0:
        fit=fit0
    elif fittype==1:
        fit=fit1
    elif fittype==2:
        fit=fit2
    elif fittype==3:
        fit=fit3
    #fit=[fit1,fit2,fit3,fit4][fittype] #doesnt work bc not all fit# defined
    end = time.time()
    Nfit=math.floor(-fit.fun)
#%% center solution (if not symmetric)
    fit_offsets=fit.x
    fit_V=constructV2(fit_offsets,fittype)
    fit_valid=np.tile(np.all(fit_V<=ewr**2,axis=0),[4,1,1]) #mask of valid points in V2
    if symmetric:
        final_offsets=fit.x
        min_diameter=2*(fit_V[fit_valid].max()**0.5+edgeexclusionwidth)
    else:
        offsets_test=fit_offsets
        def Rmax(offsets,ft,validmask):
            if fittype==0:
                V=constructV2(offsets,ft)
            else:
                offsets_test[0]=offsets #just 1 offset
                V=constructV2(offsets_test,ft)
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
        elif fittype==3:
            fit2 = optimize.minimize_scalar(Rmax,bounds=(0, (max((width,height))+yspacing)/2),args=(fittype,fit_valid))
            centered_offsets=fit_offsets
            centered_offsets[0]=fit2.x
        print(fit2)
        final_offsets=centered_offsets
        min_diameter=2*(fit2.fun**0.5+edgeexclusionwidth)
    print('\nFinal Offsets:')
    print(final_offsets)
    print('Solution diameter:')
    print(min_diameter)
    area_utilization = Nfit*width*height*4*waferdiameter**-2/np.pi
    extraoutput.set('area utilization = {:.2%}\nminimum wafer diameter = {:.1f}\nsearch time: {:.2f} sec'\
                    .format(area_utilization,math.ceil(10*min_diameter)/10,end-start))
#%% plot result
    Xcorner,Ycorner, Xextent, Yextent = CalculatePositions(final_offsets,fittype)
    #Xcen,Ycen=calcXYcen(final_offsets,fittype)
    V2=constructV2(final_offsets,fittype)
    valid=np.all(V2**0.5<=ewr,axis=0)
    partial=np.logical_and(np.any(V2**0.5<=ewr,axis=0),np.logical_not(valid))
    
    ax.clear()
    title_text='die number = {:n}   Fit Type = {}{}\nwidth = {:n}   height = {:n}\nx spacing = {:n}   y spacing = {:n}\nwafer diameter = {:n}   edge exclusion = {:n}'\
        .format(Nfit,['Uniform Grid','Shift Rows','Shift Columns','Shifted & Rotated'][fittype],['',' + Symmetric'][symmetric],width,height,xspacing,yspacing,waferdiameter,edgeexclusionwidth)
    ax.set_title(title_text,loc='left')
    
    for idx,_ in np.ndenumerate(Xcorner):
        if partial[idx]:
            ax.add_patch(Rectangle((Xcorner[idx], Ycorner[idx]),Xextent[idx], Yextent[idx],
                 edgecolor = 'xkcd:tomato', lw=1,fill=0))
    for idx,_ in np.ndenumerate(Xcorner):
        if valid[idx]:
            ax.add_patch(Rectangle((Xcorner[idx], Ycorner[idx]),Xextent[idx], Yextent[idx],
                 edgecolor = 'xkcd:blue', lw=1,facecolor='none'))#face='xkcd:light blue'      
    ax.add_patch(Circle((0,0),waferdiameter/2,edgecolor = 'k',fill=0))
    ax.add_patch(Circle((0,0),waferdiameter/2-edgeexclusionwidth,edgecolor = 'k',ls='--',fill=0))
    ax.set_aspect('equal')
    ax.autoscale()
    canvas.draw()
    button_execute.config(relief='raised')

#%% big button and outputs
button_execute=tk.Button(frame1,text='Find Best Fit',command=execute)
r+=1
button_execute.grid(column = 1, row = r,pady=30,ipadx=20,ipady=10, columnspan=2)

extraoutput= tk.StringVar()
label_extraoutput = tk.Label(frame1,textvariable = extraoutput)
r+=1
label_extraoutput.grid(column = 1, row = r, columnspan=2)

label_note = tk.Label(frame1,text = note)
r+=1
label_note.grid(column = 1, row = r,columnspan=2,rowspan=7,pady=20)

#%% excute window
root.mainloop()