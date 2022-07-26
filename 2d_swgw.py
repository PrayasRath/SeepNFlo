# -*- coding: utf-8 -*-
"""
Created on Sat May 28 10:35:14 2022

@author: praya
"""

import os
import numpy as np
import pandas as pd
import flopy
from flopy.utils import binaryfile as bf
import flopy.utils.postprocessing as pp
import matplotlib.pyplot as plt
from flopy.plot import PlotCrossSection as pcs
import math
#%%
def make_zbot(dem,aq_thick,ztype='top_thick',nlay=10,min_thick=10):
    
    nrow,ncol = dem[None,:].shape
    zbot=np.zeros((nlay,nrow,ncol)) 
    if ztype in ['top','top_thick']:
        # Original zbot approach, layers below valley bottom
        thickness=(aq_thick-min_thick)/nlay
        zbot[0]=min(dem)-min_thick
        for i in range (1,nlay):
            zbot[i,:,:]=zbot[0]-thickness*i
            
    elif ztype in ['bot','bot_thick']:
        # equal layer thickness from land surface, deepest layer thick
        
        # multiple layers above b
        lay_thick = aq_thick/nlay
         
        for i in range(nlay):
            zbot[i,:,:]=dem-(lay_thick*(i+1))
        
        zbot[-1] = zbot[-1,0,0] # flat bottom set by valley elevation and aq_thick
    elif ztype in ['spread']:
        # each column has layer thickness set by topography
        min_z = dem[0]-aq_thick
        aq_thick_array = (dem-min_z)/nlay # thickness per column
        zbot_depth = np.arange(1,nlay+1)[:,None] * aq_thick_array[None,:]    
        zbot = dem-zbot_depth
        zbot = zbot[:,None,:]
    return zbot
#%%
def dam_model(slope=None, ht=None, bottomwidth=None,topwidth=None,uhead=None,dhead=None,
                base=None,aq_thick=None,wdir=r'E:\Khosla',
                K=1.0, r_over_K=0.,ss=1E-5,sy=0.15,ncol=200, Lx=200, Ly=1e0,
                 nlay=50,nrow=1,z0=0,modelname='2d_dam', rerun=False,
                 drn_bool=None,min_thick=10,ztype='bot_thick',use_drn=False,Dam_kmult=1E-4):
    alpha=math.atan(slope)
    dx=Lx/ncol
    x=np.arange(0,Lx,dx)
    ht=ht
    top_width=topwidth/2.
    bot_width=bottomwidth/2
    basewidth=base/2
    twidth_h=int(top_width/(dx))
    bwidth_h=int(bot_width/dx)
    base_dx=int(base/2/dx)
    uh=uhead
    dh=dhead
    dam_loc = int(Lx/dx/2.)-1
    # lake_loc=int(res_loc/dx)-1 #domain world
    dem1=np.zeros(x.shape[0],dtype=np.float64)
    
    # Make reservoir depression
   # dem1[res_loc:res_loc+width_h+1] = np.linspace(z0-depth,z0,width_h+1)
    
    # Make upslope topography
    #dem1[res_loc+width_h+1:] = (x[res_loc+width_h+1:]-x[res_loc+width_h])*np.tan(alpha) + z0
    #dem1[:res_loc-width_h]=
    dem1=x*np.tan(alpha)+z0
    dem1[dam_loc-twidth_h:dam_loc+twidth_h+1]=ht
    dem1[dam_loc-bwidth_h:dam_loc-twidth_h]=np.linspace(dem1[dam_loc-bwidth_h],dem1[dam_loc-twidth_h],num=bwidth_h-twidth_h)
    dem1[dam_loc+twidth_h:dam_loc+bwidth_h]=np.linspace(dem1[dam_loc+twidth_h],dem1[dam_loc+bwidth_h],num=bwidth_h-twidth_h)
    # y=np.arange(0,(-res_loc-width_h)*dx+Lx,dx)
    # dem1[res_loc+width_h:]=y*np.tan(alpha)+dem1[res_loc+width_h-1]

    dem1=np.float32(dem1)
    plt.plot(dem1)
    
    delc=Ly/nrow
    delr=Lx/ncol
    
    cell_types = np.ones_like(x) # 1 = normal active cell
    cell_types[dam_loc-base_dx:dam_loc+base_dx+1] = 2 # 2 = dam_structure
    cell_types[0:dam_loc-base_dx] = -2 # -2 = up_stream waterbody
    cell_types[dam_loc+base_dx+1:ncol]=-3 # -3 = _stream waterbody
    
    
    cbc_fname = os.path.join(wdir,modelname+'.cbc')        
    if not os.path.isfile(cbc_fname) or rerun: # only run model if it doesn't already have outputs
        mf = flopy.modflow.Modflow(modelname, exe_name='MODFLOW-NWT_64',version='mfnwt',model_ws=wdir)
    
        
        zbot = make_zbot(dem1,aq_thick,ztype=ztype,nlay=nlay,min_thick=min_thick)
  
        nper=1
        perlen=1
        steady=[True]
        nstp=1
        # Use default for units: meters and days
        dis=flopy.modflow.ModflowDis(mf, nlay, nrow, ncol, delr=delr, 
                                     delc=delc,top=dem1[None,:],nper=nper,perlen=perlen,
                                     nstp=nstp, botm=zbot, steady=steady)

        # BAS package        
        ibound=np.ones((nlay,nrow, ncol))
#        ibound[0,:,cell_types==2]=1
        ibound[0,0,cell_types==-2]=-1 # waterbody as constant head
        ibound[0,0,cell_types==-3]=-1
        ihead =0.5*(uhead+dhead)*  np.ones((nlay, nrow, ncol), np.float)
        ihead[0,0,cell_types==-2]=uhead # waterbody constant head value
        ihead[0,0,cell_types==-3]=dhead
        bas= flopy.modflow.ModflowBas(mf, ibound=ibound, strt=ihead)
        K_arr=K*np.ones((nlay,nrow,ncol),np.float)
        K_arr[0,:,cell_types==2]=Dam_kmult
        upw = flopy.modflow.ModflowUpw(mf,laytyp=1,hk=K_arr,ss=ss,sy=sy,vka=K_arr,ipakcb=53)
        
        
#         drain=np.where(cell_types!=-2)[0]  # everywhere but the constant heads
    
#         dr_col=drain.astype(int)
#         dr_row=[0]*len(dr_col)
#         dr_layer=[0]*len(dr_col)
#         # z_down=zbot[0][0][0]
#         # dz=dem1[reservoir_bound==False]-z_down
# #        cond3=[(1e3*K*delr*delc)/(thickness/2)]*len(drx)
#         top_lay_thick = dem1[dr_col]-zbot[0,0,dr_col]
#         cond3=(K*delr*delc)/(top_lay_thick/2) # vertical conductance from cell center to drain
#         elev_vals=dem1[drain].copy()
#         idrn=np.column_stack([dr_layer,dr_row,dr_col,elev_vals,cond3])
        
#         #recharge per year
#         rech=K*r_over_K # defined by ratio
#         rech_array = rech*np.ones_like(cell_types,dtype=float)
        
#         if drn_bool is not None:
#             rech_array[drn_bool | (cell_types==-2)] = 0 # turn off recharge for cells with active drainage or constant heads 
        
#         if use_drn:
            
#             drn=flopy.modflow.ModflowDrn(mf,stress_period_data=idrn,ipakcb=53,options=['NOPRINT']) 
#             rch = flopy.modflow.ModflowRch(mf,nrchop=3,rech=rech_array[None,:],ipakcb=53)
#         else:
#             iuzfbnd = np.ones((nrow,ncol),dtype=int)
#             iuzfbnd[0,cell_types==-2] = 0
#             K_array = K*np.ones((nrow,ncol),dtype=float)
#             K_array[:,cell_types==2] = K*kmult # outside of valley
#             uzf = flopy.modflow.ModflowUzf1(mf,nuztop=1,surfdep=0.01*dx,iuzfopt=1,
#                                             vks=K_array,finf=rech_array[None,:],iuzfbnd=iuzfbnd,ipakcb=53,
#                                             )
    
        spd = {(0, 0): [ 'save head', 'save budget']}
        oc = flopy.modflow.ModflowOc(mf, stress_period_data=spd, compact=True)
        # Add nwt package to the MODFLOW model
        # pcg= flopy.modflow.ModflowPcg(mf)
        nwt = flopy.modflow.ModflowNwt(mf,headtol=0.0001,fluxtol=500,options='COMPLEX',maxiterout=2000,iprnwt=1)
    
        # Add Modpath to the model
        
        # Write the MODFLOW model input files
        mf.write_input()
    
        # Run the MODFLOW model
        success, buff = mf.run_model(silent=True,report=True)
        if buff[-1]=='  Normal termination of simulation':
                test_result='pass'
        else:
            test_result='failed'
        
        return test_result,dem1,cell_types
    else:
        return 'exists',dem1,cell_types
    
#%%
def load_wt(modelname=None,wdir=None,totim=1,
                base=None):
    hds = bf.HeadFile(os.path.join(wdir,'{}.hds'.format(modelname)))
    cbc = bf.CellBudgetFile(os.path.join(wdir,'{}.cbc'.format(modelname)))
    flf=cbc.get_data(text='FLOW LOWER FACE')[0]
    frf=cbc.get_data(text='FLOW RIGHT FACE')[0]
    
    head = hds.get_data(totim=totim)
#    qx, qy, qz = pp.get_specific_discharge((frf,None, flf), obj)
    water_table=pp.get_water_table(head,-1e+30)
#    mf = flopy.modflow.Modflow(modelname, exe_name='MODFLOW-NWT_64',version='mfnwt',model_ws=wdir)
#    gradient=pp.get_gradients(head,obj,-1e+30)
    hds.close()
    hds = None
    return water_table,head,flf,frf
# def seep_properties(modelname=None,wdir=None,
#                          cell_types=None,use_drn=False,totim=None,
#                          dem=None,dx=None,plot_bool=False,
#                          seep_depth_threshold=-1E-2)):
    
def plot_model(modelname=None,wdir=None,ax=None):
    
    mf = flopy.modflow.Modflow.load('{}.nam'.format(modelname),model_ws=wdir,exe_name='mfnwt')
    xsec = pcs(mf,line={'row':0},ax=ax)
    xsec.plot_grid(colors='k')
    # xsec.plot_bc('DRN',color='r')
    
    wt,h,flf,frf=load_wt(modelname,wdir)
    csa = xsec.plot_array(h,head=h,masked_values=[-1e+30])
    xsec.plot_surface(wt[None,:],color='b')

    xsec.plot_ibound()
    cb = plt.colorbar(csa, shrink=0.75)