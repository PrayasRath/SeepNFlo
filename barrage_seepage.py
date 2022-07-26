# -*- coding: utf-8 -*-
"""
Created on Sat May 28 10:40:29 2022

@author: praya
"""

import sys
import os
import shutil
import numpy as np
from subprocess import check_output
import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf
# Import flopy
import flopy
#%%
wdir=r'O:\Khosla\khosla_impl'
 
def Khosla_model(LX=200,LY=25,NROW=25,NCOL=200,HK=1,VKA=1,
                 d=5,HANI=1,NLAY=1,uhead=76.5,dhead=36.5,modelname='exam'):
    modelname = modelname
    ZTOP = 1.  # the "thickness" of the profile will be 1 m (= ZTOP - ZBOT)
    ZBOT = 0.
    # NLAY = 1
    # NROW = 25
    # NCOL = 200
    DELR = LX / NCOL  # recall that MODFLOW convention is DELR is along a row, thus has items = NCOL; see page XXX in AW&H (2015)
    DELC = LY / NROW  # recall that MODFLOW convention is DELC is along a column, thus has items = NROW; see page XXX in AW&H (2015)
    DELV = (ZTOP - ZBOT) / NLAY
    BOTM = np.linspace(ZTOP, ZBOT, NLAY + 1)
    
    MF = flopy.modflow.Modflow(modelname, exe_name='MODFLOW-NWT_64',version='mfnwt',model_ws=wdir)
    TOP = np.ones((NROW, NCOL),dtype=np.float)
    DIS_PACKAGE = flopy.modflow.ModflowDis(MF, NLAY, NROW, NCOL, delr=DELR, delc=DELC,
                                   top=TOP, botm=BOTM[1:], laycbd=0)
    # Variables for the BAS package
    IBOUND = np.ones((NLAY, NROW, NCOL), dtype=np.int32)  # all nodes are active (IBOUND = 1)
    
    # make the top of the profile specified head by setting the IBOUND = -1
    IBOUND[:, 0, 0:81] = -1  #don't forget arrays are zero-based!
    IBOUND[:, 0, 121:200] =-1  #don't forget arrays are zero-based!
    #IBOUND[:, 0:d, 100] = -1
    print (IBOUND)
    STRT = 40 * np.ones((NLAY, NROW, NCOL), dtype=np.float32)  # set starting head to 40 through out model domain
    STRT[:, 0, 0:81] = uhead      # pool elevation upstream of the dam
    STRT[:, 0, 121:200] = dhead     # pool elevation downstream of the dam 
    print( STRT)
    BAS_PACKAGE = flopy.modflow.ModflowBas(MF, ibound=IBOUND, strt=STRT)
    LPF_PACKAGE = flopy.modflow.ModflowLpf(MF, hk=HK, vka=VKA,chani=0,hani=HANI,ipakcb=53)
    spd = {(0, 0): [ 'save head', 'save budget']}
    OC_PACKAGE = flopy.modflow.ModflowOc(MF,stress_period_data=spd,compact=True)
    PCG_PACKAGE = flopy.modflow.ModflowPcg(MF)
    MF.write_input()
    
    silent = False  #Print model output to screen?
    pause = False   #Require user to hit enter? Doesn't mean much in Ipython notebook
    report = True   #Store the output from the model in buff
    success, buff = MF.run_model(silent=silent, pause=pause, report=report)



#Create the headfile object and grab the results for last time.
    headfile = os.path.join(wdir, modelname + '.hds')
    headfileobj = bf.HeadFile(headfile)
    cbcfile=os.path.join(wdir,modelname + '.cbc')
    cbcobj=bf.CellBudgetFile(cbcfile)
    #Get a list of times that are contained in the model
    times = headfileobj.get_times()
    print ('Headfile (' + modelname + '.hds' + ') contains the following list of times: ', times)
    frf = cbcobj.get_data(text='FLOW RIGHT FACE', totim=1)[0]
    fff = cbcobj.get_data(text='FLOW FRONT FACE', totim=1)[0]
    #FIG = plt.figure(figsize=(15,15))
    #FIG2=  plt.figure(figsize=(15,15))
    #setup contour levels and plot extent
    LEVELS = np.arange(36., 77., 1.)
    EXTENT = (DELR/2., LX - DELR/2., DELC/2., LY - DELC/2.)
    print( 'Contour Levels: ', LEVELS)
    print( 'Extent of domain: ', EXTENT)
    HEAD = headfileobj.get_data(totim=1.0)
    # #Make a contour plot on the first axis
    # AX1 = FIG.add_subplot(1, 2, 1, aspect='equal')
    # AX1.contour(np.flipud(HEAD[0, :, :]), levels=LEVELS, extent=EXTENT)
    
    # # #Make a color flood on the second axis
    # AX2 = FIG.add_subplot(1, 2, 2, aspect='equal')
    # cax = AX2.imshow(HEAD[0, :, :], extent=EXTENT, interpolation='nearest')
    # cbar = FIG.colorbar(cax, orientation='vertical', shrink=0.25)
    # # #%%
    # AX3 = FIG2.add_subplot()
    # AX3.set_xlim(75,125)
    # modelmap = flopy.plot.ModelMap(model=MF)
    # lc = modelmap.plot_grid(linewidth=1, color='gray',alpha=0.03)
    # vectors = modelmap.plot_discharge(frf, fff, head=HEAD,istep=2,jstep=2,normalize=False,color='cyan')
    # AX3.set_xlim(75,125)
    return -(min(fff[0][0])),HEAD,frf,fff,DELR,DELC,LX,LY, HEAD[0][1][121]-HEAD[0][0][121]
#%%diffrent depths
Depth=np.arange(5,505,5)
Exit_grad=[]
Head=[]
delr=[]
delc=[]
Lx=[]
DischargeY=[]
for i in range(len(Depth)):
    MF,h,flf,fff,dr,dc,x,y,eg=Khosla_model(LY=Depth[i],NROW=Depth[i])
    Exit_grad.append(eg)
    Head.append(h)
    delr.append(dr)
    delc.append(dc)
    Lx.append(x)
    DischargeY.append(MF)
#%%    
FIG = plt.figure(figsize=(15,15))
AX3 = FIG.add_subplot()
AX3.plot(np.array(Depth)/200,Exit_grad,"X-",color='k')
AX3.set_xlabel('D/L')
AX3.set_ylabel('Upward Exit Gradient [-]')
#%%
FIG = plt.figure(figsize=(15,15))
AX3 = FIG.add_subplot()
AX3.plot(np.array(Depth)/200,DischargeY,"X-",color='k')
AX3.set_xlabel('D/L')
AX3.set_ylabel('Discharge')
#%%
FIG = plt.figure(figsize=(15,15))
AX3 = FIG.add_subplot()
AX3.plot(Depth,Exit_grad,"X-")
AX3.set_xlabel('Depth(meters)')
AX3.set_ylabel('Upward Exit Gradient [-]')
LEVELS = np.arange(36., 77., 1.)
Fig3=plt.figure(figsize=(15,15))
d=[0,4,9,19]
#
EXTENT = (delr[0]/2, Lx[0] - delr[0]/2, delc[0]/2, Depth[0] - delc[0]/2)
AX1 = Fig3.add_subplot(2, 2, 1, aspect='equal')
AX1.contour(np.flipud(Head[0][0, :, :]), levels=LEVELS, extent=EXTENT)
AX1.set_title("D/L=0.025")
AX1.set_xlabel('Lx (meters)')
AX1.set_ylabel('Ly (meters)')
EXTENT = (delr[4]/2, Lx[4] - delr[4]/2, delc[4]/2, Depth[4] - delc[4]/2)
AX2 = Fig3.add_subplot(2, 2, 2, aspect='equal')
AX2.contour(np.flipud(Head[4][0, :, :]), levels=LEVELS, extent=EXTENT)
AX2.set_title("D/L=0.125")
AX2.set_xlabel('Lx (meters)')
AX2.set_ylabel('Ly (meters)')
EXTENT = (delr[19]/2, Lx[19] - delr[19]/2, delc[19]/2, Depth[19] - delc[19]/2)
AX3 = Fig3.add_subplot(2, 2, 3, aspect='equal')
AX3.contour(np.flipud(Head[9][0, :, :]), levels=LEVELS, extent=EXTENT)
AX3.set_title("D/L=0.5")
AX3.set_xlabel('Lx (meters)')
AX3.set_ylabel('Ly (meters)')
EXTENT = (delr[99]/2, Lx[99] - delr[99]/2, delc[99]/2, Depth[99] - delc[99]/2)
AX4 = Fig3.add_subplot(2, 2, 4, )
AX4.contour(np.flipud(Head[19][0, :, :]), levels=LEVELS, extent=EXTENT)
AX4.set_title("D/L=2.5")
AX4.set_xlabel('Lx (meters)')
AX4.set_ylabel('Ly (meters)')
Fig3.subplots_adjust(right=0.9)
#cbar=fig.colorbar(im,ax=[ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9],orientation='horizontal')
import matplotlib as mpl
cmap=mpl.cm.viridis
norm= mpl.colors.Normalize(vmin=36.5, vmax=76.5)
cbar_ax = Fig3.add_axes([0.92, 0.15, 0.025, 0.7])
cbar=mpl.colorbar.ColorbarBase(cbar_ax,cmap=cmap,norm=norm)#orientation='horizontal')
#cbar=fig.colorbar(im,orientation='horizontal')
#cbar.set_ticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
cbar.set_label("Head [meter]")
plt.show()
#%%
log_ani=np.arange(-5,5.5,0.1,dtype=float)
Exit_grad=[]
Head=[]
delr=[]
delc=[]
Lx=[]
Ly=[]
DischargeY=[]
for i in range(len(log_ani)):
    MF,h,flf,fff,dr,dc,x,y,eg=Khosla_model(HANI=10**(log_ani[i]))
    Exit_grad.append(eg)
    Head.append(h)
    delr.append(dr)
    delc.append(dc)
    Lx.append(x)
    Ly.append(y)
    DischargeY.append(MF)
#%% 
Fig3=plt.figure(figsize=(15,15))   
EXTENT = (delr[30]/2, Lx[30] - delr[30]/2, delc[30]/2, Ly[30] - delc[30]/2)
AX1 = Fig3.add_subplot(2, 2, 1, aspect='equal')
AX1.contour(np.flipud(Head[30][0, :, :]), levels=LEVELS, extent=EXTENT)
AX1.set_title("Ky/Kx=0.01")
AX1.set_xlabel('Lx (meters)')
AX1.set_ylabel('Ly (meters)')
EXTENT = (delr[40]/2, Lx[40] - delr[40]/2, delc[40]/2, Ly[40] - delc[40]/2)
AX2 = Fig3.add_subplot(2, 2, 2, aspect='equal')
AX2.contour(np.flipud(Head[40][0, :, :]), levels=LEVELS, extent=EXTENT)
AX2.set_title("Ky/Kx=0.1")
AX2.set_xlabel('Lx (meters)')
AX2.set_ylabel('Ly (meters)')
EXTENT = (delr[60]/2, Lx[60] - delr[60]/2, delc[60]/2, Ly[60] - delc[60]/2)
AX3 = Fig3.add_subplot(2, 2, 3, aspect='equal')
AX3.contour(np.flipud(Head[60][0, :, :]), levels=LEVELS, extent=EXTENT)
AX3.set_title("Ky/Kx=10")
AX3.set_xlabel('Lx (meters)')
AX3.set_ylabel('Ly (meters)')
EXTENT = (delr[70]/2, Lx[70] - delr[70]/2, delc[70]/2, Ly[70] - delc[70]/2)
AX4 = Fig3.add_subplot(2, 2, 4, aspect='equal')
AX4.contour(np.flipud(Head[70][0, :, :]), levels=LEVELS, extent=EXTENT)
AX4.set_title("Ky/Kx=100")
AX4.set_xlabel('Lx (meters)')
AX4.set_ylabel('Ly (meters)')
Fig3.subplots_adjust(right=0.9)
#cbar=fig.colorbar(im,ax=[ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9],orientation='horizontal')
import matplotlib as mpl
cmap=mpl.cm.viridis
norm= mpl.colors.Normalize(vmin=36.5, vmax=76.5)
cbar_ax = Fig3.add_axes([0.92, 0.15, 0.025, 0.7])
cbar=mpl.colorbar.ColorbarBase(cbar_ax,cmap=cmap,norm=norm)#orientation='horizontal')
#cbar=fig.colorbar(im,orientation='horizontal')
#cbar.set_ticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
cbar.set_label("Head [meter]")
plt.show()
#%%
FIG = plt.figure(figsize=(15,15))
AX3 = FIG.add_subplot()
AX3.plot(log_ani,Exit_grad,"X-",color='k')
AX3.set_xlabel('log10(Ky/Kx)')
AX3.set_ylabel('Upward Exit Gradient')
#%%
FIG = plt.figure(figsize=(15,15))
AX3 = FIG.add_subplot()
AX3.plot(log_ani,DischargeY,"X-",color='k')
AX3.set_xlabel('log10(Ky/Kx)')
AX3.set_ylabel('Discharge (cu.meters/second)')