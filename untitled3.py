# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 21:23:56 2022

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
 
def Khosla_modelWD(LX=200,LY=250,HK=1.0,VKA=1.0,dam_l=40.0,
                 d=0.0,HANI=1.0,NLAY=1,uhead=76.5,dhead=36.5,modelname='exam3',wdir=r'O:\Khosla\khosla_impl'):
    modelname = modelname
    ZTOP = 1.  # the "thickness" of the profile will be 1 m (= ZTOP - ZBOT)
    ZBOT = 0.
    loc=LX/2
    # NLAY = 1
    NROW = LX
    NCOL = LY
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
    IBOUND[:, 0, 0:int(((LX/2-dam_l/2)/DELR)+1)] = -1  #don't forget arrays are zero-based!
    IBOUND[:, 0, int(((LX/2+dam_l/2)/DELR)+1):NCOL] =-1  #don't forget arrays are zero-based!
    IBOUND[:, 0:int(d/DELC),int( loc/DELR)] = 0
    print (IBOUND)
    STRT = 40 * np.ones((NLAY, NROW, NCOL), dtype=np.float32)  # set starting head to 40 through out model domain
    STRT[:, 0, 0:int(((LX/2-dam_l/2)/DELR)+1)] = uhead      # pool elevation upstream of the dam
    STRT[:, 0, int(((LX/2+dam_l/2+1)/DELR)):NCOL] = dhead     # pool elevation downstream of the dam 
    #STRT[:, 0:d, 100] =40
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
    FIG = plt.figure(figsize=(15,15))
    FIG2=  plt.figure(figsize=(15,15))
    #setup contour levels and plot extent
    LEVELS = np.arange(36., 77., 1.)
    EXTENT = (DELR/2., LX - DELR/2., DELC/2., LY - DELC/2.)
    print( 'Contour Levels: ', LEVELS)
    print( 'Extent of domain: ', EXTENT)
    HEAD = headfileobj.get_data(totim=1.0)
    #Make a contour plot on the first axis
    AX1 = FIG.add_subplot(1, 2, 1, aspect='equal')
    AX1.contour(np.flipud(HEAD[0, :, :]), levels=LEVELS, extent=EXTENT)
    
    # # #Make a color flood on the second axis
    AX2 = FIG.add_subplot(1, 2, 2, aspect='equal')
    cax = AX2.imshow(HEAD[0, :, :], extent=EXTENT,vmin=dhead,vmax=uhead, interpolation='nearest')
    cbar = FIG.colorbar(cax, orientation='vertical', shrink=0.25)
    # #%%
    AX3 = FIG2.add_subplot()
    AX3.set_xlim(75,125)
    modelmap = flopy.plot.PlotMapView(model=MF)
    lc = modelmap.plot_grid(linewidth=1, color='gray',alpha=0.03)
    vectors = modelmap.plot_vector(frf, -fff,istep=10,jstep=2,normalize=True,color='cyan')
#    AX3.set_xlim(75,125)
    return -(min(fff[0][0])), (HEAD[0][1][int((loc+dam_l/2)/DELR+1)]-HEAD[0][0][int(((loc+dam_l/2)/DELR)+1)])/DELC