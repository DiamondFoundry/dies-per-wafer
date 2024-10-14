import math
import tkinter as tk

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.patches import Circle, Rectangle

from dies_per_wafer_calculator import DiesPerWaferCalculator

plt.ioff()

def run_ui():  
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

        dpw = DiesPerWaferCalculator(
            width=width,
            height=height,
            xspacing=xspacing,
            yspacing=yspacing,
            waferdiameter=waferdiameter,
            edgeexclusionwidth=edgeexclusionwidth,
            ft_grid=ft_grid,
            searchdepth=searchdepth,
            symmetric=symmetric,
            ft_ShiftRows=ft_ShiftRows, 
            ft_ShiftCols=ft_ShiftCols,
            ft_ShiftRot=ft_ShiftRot
        )

        dpw.fit()

        area_utilization = dpw.Nfit*width*height*4*waferdiameter**-2/np.pi
        extraoutput.set('area utilization = {:.2%}\nminimum wafer diameter = {:.1f}\nsearch time: {:.2f} sec'\
            .format(area_utilization,math.ceil(10*dpw.final_diameter)/10,dpw.end-dpw.start))

        if all(~np.array([ft_grid,ft_ShiftRows,ft_ShiftCols,ft_ShiftRot])):
            ax.clear()
            canvas.draw()
            extraoutput.set('No Fit Types Selected')
            button_execute.config(relief='raised')
            return
        

        valid=np.all(dpw.V2**0.5<=dpw.ewr,axis=0)
        partial=np.logical_and(np.any(dpw.V2**0.5<=dpw.ewr,axis=0),np.logical_not(valid))
        
        ax.clear()
        title_text='die number = {:n}   Fit Type = {}{}\nwidth = {:n}   height = {:n}\nx spacing = {:n}   y spacing = {:n}\nwafer diameter = {:n}   edge exclusion = {:n}'\
            .format(dpw.Nfit,['Uniform Grid','Shift Rows','Shift Columns','Shifted & Rotated'][dpw.fittype],['',' + Symmetric'][symmetric],width,height,xspacing,yspacing,waferdiameter,edgeexclusionwidth)
        ax.set_title(title_text,loc='left')

        for idx,_ in np.ndenumerate(dpw.Xcorner):
            if partial[idx]:
                ax.add_patch(Rectangle((dpw.Xcorner[idx], dpw.Ycorner[idx]),dpw.Xextent[idx], dpw.Yextent[idx],
                    edgecolor = 'xkcd:tomato', lw=1,fill=0))
            if valid[idx]:
                ax.add_patch(Rectangle((dpw.Xcorner[idx], dpw.Ycorner[idx]),dpw.Xextent[idx], dpw.Yextent[idx],
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


if __name__ == "__main__":
    run_ui()
