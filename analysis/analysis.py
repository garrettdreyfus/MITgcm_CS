from matlabglib import GLIBfromFile
from sklearn.linear_model import LinearRegression
import numpy as np
from jmd95 import dens
from jmd95_coef import rho_s_t
from scipy.integrate import quad
import gsw
import matplotlib.pyplot as plt
from datainput import *
import matplotlib
from xmitgcm import open_mdsdataset
from scipy.stats import pearsonr
rng = np.random.default_rng(0)
from tqdm import tqdm
import os
import gsw
from sympy import Symbol
from sympy import solve,nsolve,re
import ipdb
import contourpy as cpy

## Get temperature at depth)
def intTemp(depth,zgl,fname,fixsal=None):
    variables = grabMatVars(fname,('tNorth','tEast','sEast','zz','pp'))
    ## forcing temperature along eastern boundary
    tEast = np.asarray(variables["tEast"])#[0]+1.8
    sEast = np.asarray(variables["sEast"])#[0]+1.8
    zz = np.asarray(variables["zz"])[0]
    pp = np.asarray(variables["pp"])[0]
    ## forcing temperature at north eastern end of domain
    tEast = tEast[int(tEast.shape[0]-1)]
    sEast = sEast[int(sEast.shape[0]-1)]

    if fixsal:
        Tf= (0.0901-0.0575*fixsal) - (7.61*10**(-4))*pp
    else:
        Tf= (0.0901-0.0575*sEast) - (7.61*10**(-4))*pp
    print(zz[0])
    tEast = tEast-Tf
    #tEast = tEast+1.9
    f_interp = lambda xx: np.interp(xx, zz[::-1], tEast[::-1])
    results = []
    ls = []
    ## integrate and average temperature 25 meters above hub depth
    result = quad(f_interp,zgl,min(zgl+100,0), points = zz[::-1])[0]
    result = ((result/min(100,abs(zgl))))
    return result

## Get temperature at depth)
def fpAtGl(zgl,salt):
    Tf= (0.0901-0.0575*salt) - (7.61*10**(-4))*abs(zgl)
    return Tf


## Calculates slope of ice shelf from either the model parameters (param option) or from a point on the ice shelf
    # the ice shelf is linear so these methods are identical
def volume(fname):
    variables = grabMatVars(fname,('icedraft',"h","YY","xx","yy","Yicefront","XX"))
    icedraft = np.asarray(variables["icedraft"])
    h = np.asarray(variables["h"])
    icedraft[icedraft==0]=np.nan
    h = np.abs(h)
    icedraft = np.abs(icedraft)
    return np.nansum((h-icedraft)[h>icedraft]),np.nanmean((h-icedraft)[h>icedraft]),np.nansum(np.nansum([h>icedraft]))



def slope(fname,method="param"):
    if method == "lstsq":
        variables = grabMatVars(fname,('icedraft',"h","YY","xx","yy","Yicefront","XX"))
        icedraft = np.asarray(variables["icedraft"]).T
        h = np.asarray(variables["h"])
        xx = np.asarray(variables["xx"])
        yy = np.asarray(variables["yy"])
        X,Y = np.meshgrid(xx,yy)

        icedraft[icedraft==0]=np.nan
        X=X[~np.isnan(icedraft)]
        Y=Y[~np.isnan(icedraft)]
        flatclipped=icedraft[~np.isnan(icedraft)]
        A = np.vstack([X,Y, np.ones(len(X))]).T
        m1,m2, c = np.linalg.lstsq(A, flatclipped, rcond=None)[0]
        m1=np.abs(m1)
        m2=np.abs(m2)
        return np.sqrt(m1**2+m2**2)
    if method == "param":
        variables = grabMatVars(fname,('Zcdw_pt_shelf','icedraft','tEast','zz','yy',"xx","Yicefront","Hicefront","Y"))
        Y = np.asarray(variables["Y"])
        yicefront = np.asarray(variables["Yicefront"])[0][0]
        H = np.asarray(variables["Hicefront"])[0][0]
        icedraft = np.asarray(variables["icedraft"])
        Y = Y.flatten()
        icedraft = icedraft.flatten()
        ygl = np.nanmean(Y[np.nanargmax(np.abs(icedraft))])
        print(ygl)
        zgl = np.nanmin(icedraft)-np.nanmax(icedraft[icedraft!=0])
        return ((abs(zgl)-H)/(yicefront-ygl))

    if method == "grad":
        variables = grabMatVars(fname,('icedraft',"h","YY","xx","yy","Yicefront"))
        icedraft = np.asarray(variables["icedraft"])
        h = np.asarray(variables["h"])
        YY = np.asarray(variables["YY"])
        yy = np.asarray(variables["yy"])[0]
        xx = np.asarray(variables["xx"])[0]
        diff = np.abs(h)-np.abs(icedraft)
        grad = np.gradient(icedraft)
        grad[0] = (grad[0]/np.gradient(xx)[10])**2
        grad[1] = (grad[1]/np.gradient(yy)[10])**2
        #grad = np.sqrt(grad[0] + grad[1])
        grad = np.sqrt(np.sum(grad,axis=0))
        return np.nanmedian(grad[np.logical_and(icedraft!=0,diff!=0)])#np.mean(diff[np.logical_and(icedraft!=0,diff!=0)]) #+ abs(zglib-zgl)/y


def FStheory(fname,xval,regime = "default",shelf=False):

    #pull in timeseries data for returning the diagnosed meltrate 
    data = timeSeries(fname)

    #Calculate HUB from model setup file
    hub = GLIBfromFile(matVarsFile(fname))

    #We care about the mean of the model output
    for k in data.keys():
        if k != "ts":
            try:
                data[k] = np.nanmean(data[k][data["ts"]*(10**9)>7])
            except:
                data[k]=np.nan

    ##Pull in relevant geometric parameters and hydrographic forcings
    variables = grabMatVars(fname,("Hshelf","Xeast","Xwest","randtopog_height","Yicefront","Zcdw_pt_South","Zcdw_pt_shelf","tEast","sEast","icedraft","zz","h","yy","xx","saltfluxvals"))

    icedraft = np.asarray(variables["icedraft"])
    h = np.asarray(variables["h"])
    hShelf = float(abs(np.asarray(variables["Hshelf"])[0][0])-abs(200))
    #
    ##Temperature and salinity at the northern boundary
    tNorth = np.asarray(variables["tEast"])[-1,:]
    sNorth = np.asarray(variables["sEast"])[-1,:]

    vol,thick, area = volume(fname)
    ## crude but accurate way to calculate the grounding line depth
    zgl = np.nanmin(icedraft)


    ## Grab temperature at HUB depth
    ## ice shelf slope
    ices = slope(fname)
    ds = open_mdsdataset(fname,prefix=["SHIfwFlx"])

    dx = np.gradient(ds.XG)[0]
    dy = np.gradient(ds.YG)[0] 

    saltfluxvals = np.asarray(variables["saltfluxvals"])

    yice = np.asarray(variables["Yicefront"])[0][0]

    meltsaltflux = np.sum(ds.SHIfwFlx,axis=[1,2]).values*dx*dy*34.5
    polynaflux = np.sum(saltfluxvals*dx*dy)

    zz = np.asarray(variables["zz"])[0]


    # f is defined in the model setup
    f = 1.3*10**-4
    rho0 = 1025
    rhoi = 910
    Cp = 4186
    If = 334000
    #gprime_ext = data["gprime"]
    zz_i = np.linspace(np.min(zz),np.max(zz),num=200)
    sNorth_i = np.interp(zz_i,zz[::-1],sNorth[::-1])
    tNorth_i = np.interp(zz_i,zz[::-1],tNorth[::-1])
    tNorth = tNorth_i
    sNorth = sNorth_i
    zz=zz_i
    rho0 = 1025
    rhoi = 910
    Cp = 4186
    If = 334000
    W0 =  100000#(rho0*Cp)/(rhoi*If*C)
    print(fname,data['overturning_connect'])
    if regime=="default":
        if data['overturning_connect']>0.015:
            regime="connected"
        else:
            regime="disconnected"
    if regime=="connected":
        ## calculation of gprime
        if shelf:
            ds = open_mdsdataset(fname,prefix=["SHIfwFlx","SALT","THETA"])
            yindex = np.argmin(np.abs(ds.YC.values-yice-2000))
            hfac = ds.hFacC[:,yindex:yindex+2,75:125].mean(axis=1).mean(axis=1)
            cutoff = np.where(hfac.values!=1)[0][0]-1
            sNorth = ds.SALT[17:].mean(axis=0)[:,yindex:yindex+2,75:125].mean(axis=1).mean(axis=1)
            tNorth = ds.THETA[17:].mean(axis=0)[:,yindex:yindex+2,75:125].mean(axis=1).mean(axis=1)

            zz = np.asarray(variables["zz"])[0]

            zz_i = np.linspace(zz[cutoff],np.max(zz),num=200)
            sNorth_i = np.interp(zz_i,zz[cutoff::-1],sNorth[cutoff::-1])
            tNorth_i = np.interp(zz_i,zz[cutoff::-1],tNorth[cutoff::-1])
            zz = zz_i

            tNorth = tNorth_i
            sNorth = sNorth_i

            # sNorth = ds.SALT.sel(Y
        localdens = dens(sNorth,tNorth,abs(zz))
        ## density gradient
        gradd = np.abs(np.diff(localdens)/np.diff(zz))
        #average depth of all above 80th percentile
        tcline_height=np.mean(zz[:-1][gradd>np.quantile(gradd,0.85)])#+75
        zpyci = np.argmin(np.abs(np.abs(zz)-abs(tcline_height)))
        localdens = dens(sNorth,tNorth,abs(zz[zpyci]))
        rho_1i = np.logical_and(zz<zz[zpyci],zz>zz[zpyci]-50)
        rho_2i = np.logical_and(zz<zz[zpyci]+50,zz>zz[zpyci])

        meltsaltflux = np.sum(ds.SHIfwFlx,axis=[1,2]).values*dx*dy*34.5

        # beta = rho_s_t(34.5,-1.9,abs(tcline_height))[0]
        # yy = np.asarray(variables["yy"])
        # Yicefront = np.asarray(variables["Yicefront"])[0][0]
        # frontindex = np.argmin(np.abs(yy-Yicefront))
        # depth = np.abs(np.asarray(variables["h"]))
        # Hents = depth[:,frontindex][depth[:,frontindex]>200]
        # drho = (rho0/9.8)*(f*(Btotal))**(1/2)/np.nanmean(Hents)

        gprime_ext = 9.8*(np.nanmean(localdens[rho_1i])-np.nanmean(localdens[rho_2i]))/np.mean(localdens[np.logical_or(rho_1i,rho_2i)])

        Cp = 4186
        If = 334000
        deltaH = -(abs(tcline_height)- abs(hub))
        Tcdw = intTemp(hub,zgl,fname)
        stats = {"thickness":deltaH,"thermal":Tcdw,"gprime":gprime_ext,"ices":ices,"regime":regime}
        # return np.nan,np.nan,stats
        print(fname,stats)
        C =  ((rho0*Cp)/(rhoi*If*W0))

        print(stats)

        return (C)*ices*gprime_ext*deltaH/f*(Tcdw)*(60*60*24*365),-data["shiflx"],stats
    else:
        ## Let us walk through this because it requires some precisionn
        #x is the area average melt rate in m/yr
        x = Symbol('x',real=True,positive=True)

        # we convert the melt rate to a buoyancy (note the conversion to m/s of melt)
        meltflux = (x/(60*60*24*365))*(920.0)
        # melt is positive, polyna flux is negative

        Btotal = ((meltflux*area*dx*dy*34.5)+np.mean(polynaflux))/rho0*9.8*(7.8*10**(-4))
 
        # defining some useful geometric paramters, the ice draft is 200m
        z_midshelf = (zgl+(abs(zgl)-200)/2)
        depth = np.abs(np.asarray(variables["h"]))
        yy = np.asarray(variables["yy"])
        Yicefront = np.asarray(variables["Yicefront"])[0][0]
        frontindex = np.argmin(np.abs(yy-Yicefront))

        Hents = depth[:,frontindex][depth[:,frontindex]>200]

        Socean = np.mean(sNorth_i[np.abs(zz)<np.abs(np.nanmean(Hents))+200])
        Tocean = -1.9#tNorth_i[0]
        beta = rho_s_t(Socean,-1.9,abs(z_midshelf))[0]

        ## The equations are in the SI but this solves for g'_dc
        rhomin,rhomax = (rho0/9.8)*(f*(-Btotal))**(1/2)/(np.nanmax(Hents)),(rho0/9.8)*(f*(-Btotal))**(1/2)/(np.nanmin(Hents))
        rhomean = (rho0/9.8)*(f*(-Btotal))**(1/2)/(np.nanmean(Hents))
        stratterm = rhomax-rhomin
        Spolyna = Socean + rhomean/beta
        Tf = fpAtGl(z_midshelf,Spolyna)
        Tpolyna = Tocean
        D = (1-(Cp/If)*(Tf-(Tpolyna))/4)
        gprime = (9.8/1027)*(Spolyna*(1-1/D)*beta + (stratterm/6))
        stats = {}

        # 0.2 is the alpha calculated from regresssing onto these same equations but setting Bpolyna to the diagnosed Bpolyna
        alpha = 0.25
        C =  (alpha*(rho0*Cp)/(rhoi*If*W0))
        # This freezing temperature is not the bulk layer temp used in the g' calculation but rather the temperature of water leaving the cavity
        Tf = fpAtGl(200,Spolyna)

        # Remember x is m/yr, however our equations for \dot{m} are in m/s so we need the 60*60*24*365 conversion factor. Remember ice draft is 200
        solveexpr = x - (C*(60*60*24*365))*gprime*ices*((Tpolyna-Tf))*np.mean(Hents-200)/f
        disconnectedmelt = float(re(nsolve(solveexpr,x,1)))

        stats = {"thickness":np.mean(Hents-200),"thermal":((Tpolyna-Tf)).subs({x:disconnectedmelt}),"gprime":gprime.subs({x:disconnectedmelt}),"ices":ices,"regime":regime}
        # stats = {"thickness":np.mean(Hents-200),"thermal":((Tpolyna-Tf)),"gprime":gprime,"ices":ices,"regime":regime}
 
        print(stats)
        return disconnectedmelt/alpha,-data["shiflx"],stats
        # return np.nan,np.nan,stats
#condstructing depth from depth differences
def depthFromdZ(ds):
    U = ds.UVEL.values[0,:,:,:]
    fZ = list(ds.Z)
    DZ = np.asarray(fZ)
    Z = np.repeat(DZ,U.shape[1]*U.shape[2])
    Z = Z.reshape(U.shape[0],U.shape[1],U.shape[2])
    z = np.concatenate((ds.hFacC.values[:-1,:,:]-ds.hFacC.values[1:,:,:],ds.hFacC[-1:,:,:]),axis=0)
    z[z<0]=0
    z[-1,:,:]=0
    zmask = z
    return np.sum(zmask*Z,axis=0)

#Slightly grotesque way to get grid cells closest to topography. Only used for graphing.
def bottomMask(fname,ds,thresh=10):
    vals = grabMatVars(fname,("h","icedraft"))
    h = np.abs(np.asarray(vals["h"]))
    Znew = np.zeros(ds.THETA.shape[1:])
    for l in range(Znew.shape[0]):
        Znew[l,:,:]=ds.Z.values[l]
    bottom_dist = Znew-(-h.T)
    return np.logical_and(np.logical_and(bottom_dist < 50,bottom_dist>0),ds.hFacC.values>0)

#Slightly grotesque way to get grid cells closest to ice. Only used for graphing.
def icemask(fname,ds,thresh=10):
    vals = grabMatVars(fname,("h","icedraft"))
    icedraft = np.abs(np.asarray(vals["icedraft"]))
    Znew = np.zeros(ds.THETA.shape[1:])
    for l in range(Znew.shape[0]):
        Znew[l,:,:]=ds.Z.values[l]
    bottom_dist = Znew-(-icedraft.T)

    return np.logical_and(np.logical_and(bottom_dist > -20,bottom_dist<0),ds.hFacC.values>0)

def mixedLayerQuant(ds,fname):
    THETA = ds.THETA.values
    SALT = ds.SALT.values
    UVEL = ds.UVEL.values
    VVEL = ds.VVEL.values
    vals = grabMatVars(fname,("icedraft","Zcdw_pt_South"))
    icedraft = np.abs(np.asarray(vals["icedraft"]))
    max_height = float(vals["Zcdw_pt_South"][0][0])
    tcline_height = (max_height-75)/2.0+75
    idmask=np.logical_and(icedraft>0,icedraft>0)
    gprimes = []
    ssurf = []
    scdw = []
    tsurf = []
    tcdw = []
    gprime_map = np.empty([THETA.shape[0],THETA.shape[2],THETA.shape[3]])
    gprime_map[:] = 0

    Zfull = np.asarray(list(ds.Z))
    for t_index in tqdm(range(THETA.shape[0])):
        gprimes_t = []
        ssurf_t = []
        scdw_t = []
        tsurf_t = []
        tcdw_t = []
        froude_map = np.full_like(THETA[0,0,:,:],np.nan)
        for y_index in range(THETA.shape[2]):
            for x_index in range(THETA.shape[3]):
                tcast = THETA[t_index,:,y_index,x_index]
                scast = SALT[t_index,:,y_index,x_index]
                ucast = UVEL[t_index,:,y_index,x_index]
                vcast = VVEL[t_index,:,y_index,x_index]
                if (np.sum(abs(scast)>0.1)>2) and idmask[x_index,y_index]:
                    t = tcast[scast>0.1]
                    Z = Zfull[scast>0.1]
                    u = ucast[scast>0.1]
                    v = vcast[scast>0.1]
                    s = scast[scast>0.1]
                    #t = tcast
                    #s = scast
                    if np.sum(abs(s-s[0])>0.05)>0:#and np.sum(t>0.5)>0:
                        mldi = np.where(abs(s-s[0])>0.05)[0][0]
                        d = dens(s,t,Z[mldi])
                        #cdwi = np.where(t>0)[0][0]
                        rho_1 = np.nanmean(d[:mldi])
                        # rho_2 = np.nanmean(d[mldi:min(mldi*2,len(d)-1)])
                        rho_2 = np.nanmean(d[mldi:])
                        # gprimes_t.append(np.nanmax(np.abs(np.diff(d)/Z[mldi])))
                        gprimes_t.append(9.8*(rho_2-rho_1)/np.mean((rho_1,rho_2)))
                        gprime_map[t_index,y_index,x_index] = abs(Z[0]-Z[mldi])#gprimes_t[-1]
                        # scdw_t.append(np.nanmean(s[mldi:max(mldi*2,len(s)-1)]))
                        scdw_t.append(np.nanmean(s[mldi:]))
                        ssurf_t.append(np.nanmean(s[:mldi]))
                        tcdw_t.append(np.nanmean(t[mldi:]))
                        tsurf_t.append(np.nanmean(t[:mldi]))
                        velmagcdw = np.nanmean(np.sqrt(u[mldi:max(mldi*2,len(t)-1)]**2+v[mldi:max(mldi*2,len(t)-1)]**2))
                        velmagsurf = np.nanmean(np.sqrt(u[:mldi]**2+v[:mldi]**2))
                        
        gprimes.append(np.nanmean(gprimes_t))
        ssurf.append(np.nanmean(ssurf_t))
        scdw.append(np.nanmean(scdw_t))
        tsurf.append(np.nanmean(tsurf_t))
        tcdw.append(np.nanmean(tcdw_t))
    # plt.imshow(np.mean(gprime_map[10:],axis=0),vmin=0,vmax=0.001)
    # plt.imshow(np.mean(gprime_map[10:],axis=0),vmin=0,vmax=100)
    # plt.colorbar()
    # plt.show()
    return gprimes,ssurf,scdw,tsurf,tcdw


def timeSeries(fname,refresh=False):
    print(fname)
    fnameparts = fname.split("/")
    slug = fnameparts[-3]+"-"+fnameparts[-2]+".pickle"

    if os.path.isfile("data/modelTimeSeries/"+slug) and not refresh:
        with open("data/modelTimeSeries/"+slug,"rb") as f:
            output = pickle.load(f)
        return output

    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,prefix=["THETA","SALT","momKE","SHIfwFlx","VVEL","UVEL","WVEL","RHOAnoma","MXLDEPTH"],ignore_unknown_vars=True,iters=times,extra_variables = extra_variables)
    ## theta plot
    ts = ds.time.values*grabDeltaT(fname)/60.0/60.0/24.0/365.0
    tsnew = np.full_like(ts,0,dtype=float)
    tsnew[:] = ts
    ts = tsnew
    ts = ts/1000000000
    times=ts
    yice = grabMatVars(fname,("Yicefront"))["Yicefront"][0][0]
    yshelf = grabMatVars(fname,("Yshelfbreak"))["Yshelfbreak"][0][0]
    ycoast = grabMatVars(fname,("Ycoast"))["Ycoast"][0][0]
    Xeast = grabMatVars(fname,("Xeast"))["Xeast"][0][0]
    Xwest = grabMatVars(fname,("Xwest"))["Xwest"][0][0]
    
    variables = grabMatVars(fname,("tEast",'sEast','zz'))
    tNorth = np.asarray(variables["tEast"])[-1,:]
    sNorth = np.asarray(variables["sEast"])[-1,:]

    cavity = ds.where(ds.YC<yice-10000,drop=True)
    cavity = cavity.where(cavity.XC<np.max(cavity.XC.values)-10000,drop=True)
    cavity = cavity.where(cavity.XC>np.min(cavity.XC.values)+10000,drop=True)
       

    volume = (ds.hFacC * ds.drF * ds.rA).values
    cavvolume = (cavity.hFacC * cavity.drF * cavity.rA).values
    volumesum = np.sum(volume)
    cavvolumesum = np.sum(cavvolume)

    thetas=[]
    THETA = cavity.THETA.values
    for k in range(cavity.THETA.shape[0]):
        thetas.append((np.sum(THETA[k]*cavvolume))/cavvolumesum)
        #thetas.append(np.mean(THETA[k][THETA[k]!=0]))

    salts = []
    SALT = cavity.SALT.values

    for k in range(cavity.THETA.shape[0]):
        salts.append((np.sum(SALT[k]*cavvolume))/cavvolumesum)

    kes = []
    momKE = cavity.momKE.values
    for k in range(cavity.THETA.shape[0]):
        kes.append((np.sum(momKE[k]*cavvolume))/cavvolumesum)


    cavity = ds.where(ds.YC<yice,drop=True)

    incavity = []
    vvel = ds.VVEL.values
    ht = vvel*(ds.THETA.values+1.8)*(ds.RHOAnoma.values+1000)
    #frontmask = ds.hFacC[:,np.nanargmin(np.abs(ds.VVEL.YG-(yice-10000))),:]==0
    index = np.argmin(np.abs(ds.VVEL.YG.values-(yice-50000)))

    polyna = ds.where(ds.YC>=yice,drop=True)
    polynasaltdiff = polyna.where(polyna.Z>-300,drop=True).mean(dim=["Z","XC","YC"]) - polyna.where(polyna.Z<-300,drop=True).mean(dim=["Z","XC","YC"])
    polynasaltdiff = polynasaltdiff.SALT.values
    polyna = polyna.where(polyna.YC<=yice+5000,drop=True)
    mxldepth = polyna.MXLDEPTH.mean(dim=["XC","YC"]).values
    ieast = np.argmin(np.abs(Xeast-ds.XC.values))
    iwest = np.argmin(np.abs(Xwest-ds.XC.values))
    mxldepth = polyna.MXLDEPTH.values[:,:,iwest:ieast]

    zfull = polyna.Z.broadcast_like(polyna.hFacC)
    zfull.values=zfull.values*polyna.hFacC.values
    mxldiff = []
    for i in range(mxldepth.shape[0]):
    	mxldiff.append(-mxldepth[i]-zfull.min(dim="Z").values[:,iwest:ieast])

    mxldepth = -np.mean(mxldepth,axis=1)
    mxldiff = np.mean(mxldiff,axis=1)

    entrance = polyna.where(polyna.Z<-250,drop=True)
    entrancefrac = entrance.hFacC.values
    entrancesalt = []
    entrancetheta = []
    for i in range(ds.THETA.shape[0]):
        entrancesalt.append(np.sum(entrance.SALT.values[i]*entrancefrac)/np.sum(entrancefrac))
        entrancetheta.append(np.sum(entrance.THETA.values[i]*entrancefrac)/np.sum(entrancefrac))

    frontmask = ds.hFacC[:,index,:]
    sliceindex=index
    for k in range(ds.THETA.shape[0]):
        vvelcross = vvel[k,:,sliceindex,:]
        htcross = ht[k,:,sliceindex,:]*np.array(ds.drF.values)[:,None]
        #htcross[frontmask]=np.nan
        #vvelcross[frontmask]=np.nan
        incavity.append(np.nansum(htcross*frontmask))

    mask = ~np.isnan(ds.SHIfwFlx.values)
    shflx = ds.SHIfwFlx.values
    shiwflxs = []
    vals = grabMatVars(fname,("h","icedraft"))
    icedraft = np.abs(np.asarray(vals["icedraft"]))
    h = np.abs(np.asarray(vals["h"]))

    icemaskm = np.logical_and(icedraft>0,icedraft<h)
    for k in range(shflx.shape[0]):
        shiwflxs.append(np.mean(shflx[k][icemaskm.T])*(60*60*24*365)*(1/920.0))

    avgbts = [] #barotropic_streamfunction_max(fname)
    nanmask = ~np.isnan(shiwflxs)
    #np.sum(shflx*mask,axis=(1,2))*(60*60*24*365)*(1/920.0)*(1/np.sum(mask,axis=(1,2)))
    bottomtemps = thetas    # print(bottomtemps.shape)
    icem = icemask(fname,ds)
    THETA = ds.THETA.values
    SALT = ds.SALT.values
    UVEL = ds.UVEL.values
    VVEL = ds.VVEL.values
    WVEL = ds.WVEL.values
    VEL = np.sqrt(UVEL**2 + VVEL**2 + WVEL**2)
    #plt.imshow(icemaskm)
    #plt.show()
    icem[:,~icemaskm.T] = 0
    icesurfacetemps = []
    icesurfacevels = []
    icesurfacesalts = []
    meltapprox=[]
    #offshorefraction = massBoxes(fname,ds=ds)

    for k in range(THETA.shape[0]):
        icesurfacetemps.append(np.nansum(THETA[k][icem])/np.sum(icem))
        icesurfacevels.append(np.nansum(VEL[k][icem])/np.sum(icem))
        icesurfacesalts.append(np.nansum(SALT[k][icem])/np.sum(icem))
        meltapprox.append(np.nansum((THETA[k]*VEL[k])[icem])/np.sum(icem))

    gprimes,ssurf,scdw,tsurf,tcdw = mixedLayerQuant(ds,fname)
    froude = np.empty(ts.shape)
    #gprimes,ssurf,scdw,tsurf,tcdw,froude = np.empty(ts.shape),np.empty(ts.shape),np.empty(ts.shape),np.empty(ts.shape),np.empty(ts.shape),np.empty(ts.shape)
    output = {"ts":np.asarray(ts)[nanmask],"theta":np.asarray(thetas)[nanmask],"salt":np.asarray(salts)[nanmask],\
        "kes":np.asarray(kes)[nanmask],"avgbts":np.asarray([]),\
        "shiflx":np.asarray(shiwflxs)[nanmask],"bottemp":np.asarray(bottomtemps)[nanmask],"icesurfacetemp":np.asarray(icesurfacetemps)[nanmask],\
        "icesurfacesalt":np.asarray(icesurfacesalts)[nanmask],\
        "icesurfacevel":np.asarray(icesurfacevels)[nanmask],"meltapprox":np.asarray(meltapprox)[nanmask],\
        "incavity":np.asarray(incavity)[nanmask],"gprime":np.asarray(gprimes)[nanmask],\
        "ssurf":np.asarray(ssurf)[nanmask],"scdw":np.asarray(scdw)[nanmask],"tsurf":np.asarray(tsurf)[nanmask],"tcdw":np.asarray(tcdw)[nanmask],\
        "polynasaltdiff":np.asarray(polynasaltdiff)[nanmask],\
        "mxldepth":np.asarray(mxldepth)[nanmask],\
        "mxldiff":np.asarray(mxldiff)[nanmask],\
        "entrancesalt":np.asarray(entrancesalt)[nanmask],\
        "entrancetheta":np.asarray(entrancetheta)[nanmask],
        "overturning_connect":np.asarray(overturning_summary(ds,yice,ycoast))[nanmask]\
    }
    with open("data/modelTimeSeries/"+slug,"wb") as f:
        pickle.dump(output,f)
    return output


def massBoxes(fname,ds=None,show=False):

    times=getIterNums(fname)

    if not ds:
        extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
        ds = open_mdsdataset(fname,prefix=["UVEL","VVEL","THETA","SALT","momKE","SHIfwFlx","WVEL","RHOAnoma"],ignore_unknown_vars=True,iters=times,extra_variables = extra_variables)
    times = np.asarray( [0]+list(ds.time.values))*(10**-9)

    yice = grabMatVars(fname,("Yicefront"))["Yicefront"][0][0]
    yshelf = grabMatVars(fname,("Yshelfbreak"))["Yshelfbreak"][0][0]
    ycoast = grabMatVars(fname,("Ycoast"))["Ycoast"][0][0]
    xeast = grabMatVars(fname,("Xeast"))["Xeast"][0][0]
    xwest = grabMatVars(fname,("Xwest"))["Xwest"][0][0]
    saltfluxvals = grabMatVars(fname,("saltfluxvals"))["saltfluxvals"].T
    vvel = ds.VVEL.values
    #vvel = ds.SALT.values*vvel#(vvel + np.roll(vvel,-1,axis=2))/2
    uvel = ds.UVEL.values
    wvel = ds.WVEL.values
    #uvel = ds.UVEL.values
    #uvel = ds.SALT.values*uvel#(uvel + np.roll(uvel,-1,axis=3))/2
    merid_salt_transport = ds.dxG.values*(ds.hFacS*ds.drF).values*vvel*1027
    zonal_salt_transport = ds.dyC.values*(ds.hFacW*ds.drF).values*uvel*1027
    top_salt_transport = ds.rA.values*wvel*1027
    #zonal_salt_transport =  ds.SALT.values
    #fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)

    #infront_polyna_y = np.argmin(np.abs(ds.YG.values-(yice+10000)))-1
    #back_polyna_y = np.argmax(saltfluxvals[:,91]!=0)+1
    back_polyna_y = np.argmax(saltfluxvals[:,91]!=0)-1
    #back_polyna_y = np.argmin(np.abs(ds.YG.values-(yice)))+1
    #infront_polyna_y = saltfluxvals.shape[0]-np.argmax(saltfluxvals[::-1,91]!=0)-1
    infront_polyna_y = saltfluxvals.shape[0]-np.argmax(saltfluxvals[::-1,91]!=0)-2

    print(ds.YC[back_polyna_y].values)
    print(ds.YC[infront_polyna_y].values)

    #west_polyna_x = np.argmax(saltfluxvals[infront_polyna_y-1,:]!=0)
    #east_polyna_x = saltfluxvals.shape[1]-np.argmax(saltfluxvals[infront_polyna_y-1,::-1]!=0)-1

    east_polyna_x = np.argmin(np.abs(ds.XC.values-(xeast)))+1
    west_polyna_x = np.argmin(np.abs(ds.XC.values-(xwest)))

    topface_i = np.argmin(np.abs(np.abs(ds.Z.values)-250))
    top_polyna_face = top_salt_transport[:,topface_i,back_polyna_y:infront_polyna_y,west_polyna_x:east_polyna_x]

    north_polyna_face = merid_salt_transport[:,topface_i:,infront_polyna_y,west_polyna_x:east_polyna_x]
    south_polyna_face = merid_salt_transport[:,topface_i:,back_polyna_y,west_polyna_x:east_polyna_x]
    west_polyna_face = zonal_salt_transport[:,topface_i:,back_polyna_y:infront_polyna_y,west_polyna_x]
    east_polyna_face = zonal_salt_transport[:,topface_i:,back_polyna_y:infront_polyna_y,east_polyna_x]

 
    dx = np.gradient(ds.XG)[0]
    dy = np.gradient(ds.YG)[0]
    total = np.sum(south_polyna_face,axis=(1,2))+np.sum(west_polyna_face,axis=(1,2)) \
            -np.sum(east_polyna_face,axis=(1,2))-np.sum(north_polyna_face,axis=(1,2)) \
            -np.sum(top_polyna_face,axis=(1,2))
    salt = ds.SALT.values 
    saltinpolyna = []
    fluxdivvol = []
    polynaflux = -np.sum(saltfluxvals[back_polyna_y:infront_polyna_y+1,west_polyna_x:east_polyna_x+1])*dx*dy
    for k in range(salt.shape[0]):
        salttslice = salt[k]
        vol = ((ds.hFacC*ds.drF)*dx*dy)
        averagesalt = ds.SALT[k]*vol
        averagesalt = np.sum(averagesalt.values[:,back_polyna_y:infront_polyna_y,west_polyna_x:east_polyna_x],axis=(0,1,2))/np.sum(vol.values[:,back_polyna_y:infront_polyna_y,west_polyna_x:east_polyna_x],axis=(0,1,2))
    
        saltinpolyna.append(averagesalt)
        W = (-west_polyna_x+east_polyna_x)*dx
        L = (-back_polyna_y+infront_polyna_y)*dy
        if k ==0:
                deltaT = float(times[0])
        else:
                deltaT = float(times[k+1]-times[k])
        fluxdivvol.append(deltaT*(total[k]+polynaflux)/(np.sum(vol.values[:,back_polyna_y:infront_polyna_y,west_polyna_x:east_polyna_x],axis=(0,1,2))))
    fluxvol = np.cumsum(fluxdivvol)
    fluxvol = fluxvol+saltinpolyna[0]

    #topavg = -np.sum(top_polyna_face,axis=(1,2))/10**8
    #northavg = -np.sum(north_polyna_face,axis=(1,2))/10**8
    #eastavg = -np.sum(east_polyna_face,axis=(1,2))/10**8
    #southavg = np.sum(south_polyna_face,axis=(1,2))/10**8
    #westavg = np.sum(west_polyna_face,axis=(1,2))/10**8
    #top_polyna_face[top_polyna_face>0]=0
    #north_polyna_face[north_polyna_face>0]=0
    #east_polyna_face[east_polyna_face>0]=0
    ##south_polyna_face[south_polyna_face<0]=0
    #west_polyna_face[west_polyna_face<0]=0

    topavg = -np.sum(top_polyna_face,axis=(1,2))/10**8
    northavg = -np.sum(north_polyna_face,axis=(1,2))/10**8
    eastavg = -np.sum(east_polyna_face,axis=(1,2))/10**8
    southavg = np.sum(south_polyna_face,axis=(1,2))/10**8
    westavg = np.sum(west_polyna_face,axis=(1,2))/10**8

    totalavg = total/10**8
    if show:
        matplotlib.rcParams.update({'font.size': 22})
        fig, (ax1) = plt.subplots(1,1,figsize=(10,12))
        ax1.axhline(y=0,linewidth=1, color='black')
        barlabels = ['top','north','east','south','west','total']
        avgs = [np.mean(topavg[10:]),np.mean(northavg[10:]),np.mean(eastavg[10:]),np.mean(southavg[10:]),np.mean(westavg[10:]),np.mean(totalavg[10:])]
        ax1.bar(barlabels,avgs,width=0.5)
        #ax1.set_ylim(-30,30)
        ax1.set_ylabel("(g/kg)*kg/s (in 10**8)")
        plt.xticks(rotation=30, ha='right')
        plt.show()
    return topavg#/(northavg)

def overturning_summary(ds,yice,ycoast):
    vvel = ds.VVEL
    fulltransport = (vvel*ds.hFacS*ds.drF*ds.dxG)
    connect_frac = []
    for i in tqdm(range(len(fulltransport))):
        transport = fulltransport[i].mean(dim="XC").where(ds.YG<175000)
        transport = transport.cumsum(dim="Z")
        streamlow = transport.where(np.abs(transport.YG-(yice+15000))<6000).min().values
        streamhigh = transport.where(np.abs(transport.YG-(yice+15000))<6000).max().values
        levels=np.linspace(streamlow+100,streamhigh-100)
        cg = cpy.contour_generator(transport.YG,transport.Z,transport.values)
        bools = []
        for l in levels:
            if abs(l)>3:
                for seg in cg.lines(l):
                    xmax = seg[:,0].max()
                    xmin = seg[:,0].min()
                    bools.append((xmax>yice+21000) and (xmin<yice-25000))
        connect_frac.append((np.sum(bools)/len(bools)))
    print(connect_frac)
    return connect_frac



def saltBoxes(fname):

    times=getIterNums(fname)

    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    ds = open_mdsdataset(fname,prefix=["UVEL","VVEL","THETA","SALT","momKE","SHIfwFlx","VVELSLT","UVELSLT","WVEL","WVELSLT","RHOAnoma"],ignore_unknown_vars=True,iters=times,extra_variables = extra_variables)
    times = np.asarray( [0]+list(ds.time.values))*(10**-9)

    yice = grabMatVars(fname,("Yicefront"))["Yicefront"][0][0]
    yshelf = grabMatVars(fname,("Yshelfbreak"))["Yshelfbreak"][0][0]
    ycoast = grabMatVars(fname,("Ycoast"))["Ycoast"][0][0]
    xeast = grabMatVars(fname,("Xeast"))["Xeast"][0][0]
    xwest = grabMatVars(fname,("Xwest"))["Xwest"][0][0]
    saltfluxvals = grabMatVars(fname,("saltfluxvals"))["saltfluxvals"].T

    vvelslt = ds.VVELSLT.values
    uvelslt = ds.UVELSLT.values
    wvelslt = ds.WVELSLT.values

    salt = ds.SALT.values

    vvel = ds.VVEL.values
    uvel = ds.UVEL.values
    wvel = ds.WVEL.values

    merid_salt_transport = ds.dxG.values*(ds.hFacS*ds.drF).values*vvelslt*1027
    zonal_salt_transport = ds.dyC.values*(ds.hFacW*ds.drF).values*uvelslt*1027
    top_salt_transport = ds.rA.values*wvelslt*1027

    #vbarsbar = ds.dxG.values*(ds.hFacS*ds.drF).values*vvel*(salt+np.roll(salt,axis=0)*1027
    vbarsbar = ds.dxG.values*(ds.hFacS*ds.drF).values*vvel*salt*1027
    ubarsbar = ds.dyC.values*(ds.hFacW*ds.drF).values*uvel*salt*1027
    wbarsbar = ds.dyC.values*(ds.hFacW*ds.drF).values*wvel*salt*1027

    back_polyna_y = np.argmax(saltfluxvals[:,91]!=0)
    infront_polyna_y = saltfluxvals.shape[0]-np.argmax(saltfluxvals[::-1,91]!=0)-1

    east_polyna_x = np.argmin(np.abs(ds.XC.values-(xeast)))+1
    west_polyna_x = np.argmin(np.abs(ds.XC.values-(xwest)))

    topface_i = np.argmin(np.abs(np.abs(ds.Z.values)-200))
    # topface_i = np.argmin(np.abs(np.abs(ds.Z.values)-00))

    top_polyna_face = top_salt_transport[:,topface_i,back_polyna_y:infront_polyna_y,west_polyna_x:east_polyna_x-1]
    north_polyna_face = merid_salt_transport[:,topface_i:,infront_polyna_y,west_polyna_x:east_polyna_x]
    south_polyna_face = merid_salt_transport[:,topface_i:,back_polyna_y,west_polyna_x:east_polyna_x]
    west_polyna_face = zonal_salt_transport[:,topface_i:,back_polyna_y:infront_polyna_y,west_polyna_x]
    east_polyna_face = zonal_salt_transport[:,topface_i:,back_polyna_y:infront_polyna_y,east_polyna_x]

    top_polyna_face_barbar = wbarsbar[:,topface_i,back_polyna_y:infront_polyna_y,west_polyna_x:east_polyna_x-1]
    north_polyna_face_barbar = vbarsbar[:,topface_i:,infront_polyna_y,west_polyna_x:east_polyna_x]
    south_polyna_face_barbar = vbarsbar[:,topface_i:,back_polyna_y,west_polyna_x:east_polyna_x]
    west_polyna_face_barbar = ubarsbar[:,topface_i:,back_polyna_y:infront_polyna_y,west_polyna_x]
    east_polyna_face_barbar = ubarsbar[:,topface_i:,back_polyna_y:infront_polyna_y,east_polyna_x]
 
    dx = np.gradient(ds.XG)[0]
    dy = np.gradient(ds.YG)[0]
    total = np.sum(south_polyna_face,axis=(1,2))+np.sum(west_polyna_face,axis=(1,2)) \
            -np.sum(east_polyna_face,axis=(1,2))-np.sum(north_polyna_face,axis=(1,2))-np.sum(top_polyna_face,axis=(1,2))

    salt = ds.SALT.values 
    saltinpolyna = []
    fluxdivvol = []
    polynaflux = -np.sum(saltfluxvals[back_polyna_y:infront_polyna_y+1,west_polyna_x:east_polyna_x+1])*dx*dy
    for k in range(salt.shape[0]):
        salttslice = salt[k]
        vol = ((ds.hFacC*ds.drF)*dx*dy)
        averagesalt = ds.SALT[k]*vol
        averagesalt = np.sum(averagesalt.values[:,back_polyna_y:infront_polyna_y,west_polyna_x:east_polyna_x],axis=(0,1,2))/np.sum(vol.values[:,back_polyna_y:infront_polyna_y,west_polyna_x:east_polyna_x],axis=(0,1,2))
    
        saltinpolyna.append(averagesalt)
        W = (-west_polyna_x+east_polyna_x)*dx
        L = (-back_polyna_y+infront_polyna_y)*dy
        if k ==0:
                deltaT = float(times[0])
        else:
                deltaT = float(times[k+1]-times[k])
        fluxdivvol.append(deltaT*(total[k]+polynaflux)/(np.sum(vol.values[:,back_polyna_y:infront_polyna_y,west_polyna_x:east_polyna_x],axis=(0,1,2))))
    fluxvol = np.cumsum(fluxdivvol)
    fluxvol = fluxvol+saltinpolyna[0]

    barlabels = ['North','East','South','West','Top']
    northavg = np.mean(-np.sum(north_polyna_face,axis=(1,2))[15:])/10**8
    eastavg = np.mean(-np.sum(east_polyna_face,axis=(1,2))[15:])/10**8
    southavg = np.mean(np.sum(south_polyna_face,axis=(1,2))[15:])/10**8
    westavg = np.mean(np.sum(west_polyna_face,axis=(1,2))[15:])/10**8
    topavg = np.mean(-np.sum(top_polyna_face,axis=(1,2))[15:])/10**8
    totalavg = np.mean(total[15:])/10**8
    matplotlib.rcParams.update({'font.size': 22})
    fig, (ax1) = plt.subplots(1,1,figsize=(14,12))
    ax1.axhline(y=0,linewidth=1, color='black')
    avgs = [northavg,eastavg,southavg,westavg,topavg]
    ax1.bar(barlabels,avgs,width=0.5,label="Mean")
    #ax1.set_ylim(-30,30)
    ax1.set_ylabel(r'$\frac{g}{kg}~\frac{kg}{s} (10^8)$')
    plt.xticks(rotation=30, ha='right')


    barlabels = ['North','East','South','West','Top']
    northavg = np.mean(-np.sum(north_polyna_face-north_polyna_face_barbar,axis=(1,2))[10:])/10**8
    eastavg = np.mean(-np.sum(east_polyna_face-east_polyna_face_barbar,axis=(1,2))[10:])/10**8
    southavg = np.mean(np.sum(south_polyna_face-south_polyna_face_barbar,axis=(1,2))[10:])/10**8
    westavg = np.mean(np.sum(west_polyna_face-west_polyna_face_barbar,axis=(1,2))[10:])/10**8
    topavg = np.mean(-np.sum(top_polyna_face-top_polyna_face_barbar,axis=(1,2))[10:])/10**8
    avgs = [northavg,eastavg,southavg,westavg,topavg]
    ax1.bar(barlabels,avgs,width=0.25,align='center',color="red",label="Eddy")
    plt.legend()
    plt.savefig("saltboxes.png")
    plt.tight_layout()

    plt.show()
 


def saltBoxesnew(fname):

    times=getIterNums(fname)

    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    ds = open_mdsdataset(fname,prefix=["UVEL","VVEL","THETA","SALT","momKE","SHIfwFlx","VVELSLT","UVELSLT","WVEL","WVELSLT","RHOAnoma"],ignore_unknown_vars=True,iters=times,extra_variables = extra_variables)
    times = np.asarray( [0]+list(ds.time.values))*(10**-9)

    yice = grabMatVars(fname,("Yicefront"))["Yicefront"][0][0]
    yshelf = grabMatVars(fname,("Yshelfbreak"))["Yshelfbreak"][0][0]
    ycoast = grabMatVars(fname,("Ycoast"))["Ycoast"][0][0]
    xeast = grabMatVars(fname,("Xeast"))["Xeast"][0][0]
    xwest = grabMatVars(fname,("Xwest"))["Xwest"][0][0]
    saltfluxvals = grabMatVars(fname,("saltfluxvals"))["saltfluxvals"].T

    vvelslt = ds.VVELSLT.values
    uvelslt = ds.UVELSLT.values
    wvelslt = ds.WVELSLT.values

    salt = ds.SALT.values

    vvel = ds.VVEL.values
    uvel = ds.UVEL.values
    wvel = ds.WVEL.values

    merid_salt_transport = ds.dxG.values*(ds.hFacS*ds.drF).values*vvelslt*1027
    zonal_salt_transport = ds.dyC.values*(ds.hFacW*ds.drF).values*uvelslt*1027
    top_salt_transport = ds.rA.values*wvelslt*1027

    #vbarsbar = ds.dxG.values*(ds.hFacS*ds.drF).values*vvel*(salt+np.roll(salt,axis=0)*1027
    vbarsbar = ds.dxG.values*(ds.hFacS*ds.drF).values*vvel*salt*1027
    ubarsbar = ds.dyC.values*(ds.hFacW*ds.drF).values*uvel*salt*1027
    wbarsbar = ds.dyC.values*(ds.hFacW*ds.drF).values*wvel*salt*1027

    back_polyna_y = np.argmax(saltfluxvals[:,91]!=0)
    infront_polyna_y = saltfluxvals.shape[0]-np.argmax(saltfluxvals[::-1,91]!=0)-1

    east_polyna_x = np.argmin(np.abs(ds.XC.values-(xeast)))+1
    west_polyna_x = np.argmin(np.abs(ds.XC.values-(xwest)))

    topface_i = np.argmin(np.abs(np.abs(ds.Z.values)-200))
    # topface_i = np.argmin(np.abs(np.abs(ds.Z.values)-00))

    top_polyna_face = top_salt_transport[:,topface_i,back_polyna_y:infront_polyna_y,west_polyna_x:east_polyna_x-1]
    north_polyna_face = merid_salt_transport[:,topface_i:,infront_polyna_y,west_polyna_x:east_polyna_x]
    south_polyna_face = merid_salt_transport[:,topface_i:,back_polyna_y,west_polyna_x:east_polyna_x]
    west_polyna_face = zonal_salt_transport[:,topface_i:,back_polyna_y:infront_polyna_y,west_polyna_x]
    east_polyna_face = zonal_salt_transport[:,topface_i:,back_polyna_y:infront_polyna_y,east_polyna_x]

    top_polyna_face_barbar = wbarsbar[:,topface_i,back_polyna_y:infront_polyna_y,west_polyna_x:east_polyna_x-1]
    north_polyna_face_barbar = vbarsbar[:,topface_i:,infront_polyna_y,west_polyna_x:east_polyna_x]
    south_polyna_face_barbar = vbarsbar[:,topface_i:,back_polyna_y,west_polyna_x:east_polyna_x]
    west_polyna_face_barbar = ubarsbar[:,topface_i:,back_polyna_y:infront_polyna_y,west_polyna_x]
    east_polyna_face_barbar = ubarsbar[:,topface_i:,back_polyna_y:infront_polyna_y,east_polyna_x]
 
    dx = np.gradient(ds.XG)[0]
    dy = np.gradient(ds.YG)[0]
    total = np.sum(south_polyna_face,axis=(1,2))+np.sum(west_polyna_face,axis=(1,2)) \
            -np.sum(east_polyna_face,axis=(1,2))-np.sum(north_polyna_face,axis=(1,2))-np.sum(top_polyna_face,axis=(1,2))

    salt = ds.SALT.values 
    saltinpolyna = []
    fluxdivvol = []
    polynaflux = -np.sum(saltfluxvals[back_polyna_y:infront_polyna_y+1,west_polyna_x:east_polyna_x+1])*dx*dy
    for k in range(salt.shape[0]):
        salttslice = salt[k]
        vol = ((ds.hFacC*ds.drF)*dx*dy)
        averagesalt = ds.SALT[k]*vol
        averagesalt = np.sum(averagesalt.values[:,back_polyna_y:infront_polyna_y,west_polyna_x:east_polyna_x],axis=(0,1,2))/np.sum(vol.values[:,back_polyna_y:infront_polyna_y,west_polyna_x:east_polyna_x],axis=(0,1,2))
    
        saltinpolyna.append(averagesalt)
        W = (-west_polyna_x+east_polyna_x)*dx
        L = (-back_polyna_y+infront_polyna_y)*dy
        if k ==0:
                deltaT = float(times[0])
        else:
                deltaT = float(times[k+1]-times[k])
        fluxdivvol.append(deltaT*(total[k]+polynaflux)/(np.sum(vol.values[:,back_polyna_y:infront_polyna_y,west_polyna_x:east_polyna_x],axis=(0,1,2))))
    fluxvol = np.cumsum(fluxdivvol)
    fluxvol = fluxvol+saltinpolyna[0]

    barlabels = ['North','East','South','West','Top']
    northavg = np.mean(-np.sum(north_polyna_face,axis=(1,2))[15:])/10**8
    eastavg = np.mean(-np.sum(east_polyna_face,axis=(1,2))[15:])/10**8
    southavg = np.mean(np.sum(south_polyna_face,axis=(1,2))[15:])/10**8
    westavg = np.mean(np.sum(west_polyna_face,axis=(1,2))[15:])/10**8
    topavg = np.mean(-np.sum(top_polyna_face,axis=(1,2))[15:])/10**8
    totalavg = np.mean(total[15:])/10**8
    matplotlib.rcParams.update({'font.size': 22})
    fig, (ax1) = plt.subplots(1,1,figsize=(14,12))
    ax1.axhline(y=0,linewidth=1, color='black')
    avgs = [northavg,eastavg,southavg,westavg,topavg]
    ax1.bar(barlabels,avgs,width=0.5,label="Mean")
    #ax1.set_ylim(-30,30)
    ax1.set_ylabel(r'$\frac{g}{kg}~\frac{kg}{s} (10^8)$')
    plt.xticks(rotation=30, ha='right')


    barlabels = ['North','East','South','West','Top']
    northavg = np.mean(-np.sum(north_polyna_face-north_polyna_face_barbar,axis=(1,2))[10:])/10**8
    eastavg = np.mean(-np.sum(east_polyna_face-east_polyna_face_barbar,axis=(1,2))[10:])/10**8
    southavg = np.mean(np.sum(south_polyna_face-south_polyna_face_barbar,axis=(1,2))[10:])/10**8
    westavg = np.mean(np.sum(west_polyna_face-west_polyna_face_barbar,axis=(1,2))[10:])/10**8
    topavg = np.mean(-np.sum(top_polyna_face-top_polyna_face_barbar,axis=(1,2))[10:])/10**8
    avgs = [northavg,eastavg,southavg,westavg,topavg]
    ax1.bar(barlabels,avgs,width=0.25,align='center',color="red",label="Eddy")
    plt.legend()
    plt.savefig("saltboxes.png")
    plt.tight_layout()

    plt.show()
 

def breakdown(fname,xval,fig,axises,title="",color="blue",marker="o",znum=0,count=0):
    ((ax1,ax2),(ax3,ax4)) = axises
    _,_,stats = FStheory(fname,None)
    if "thickness" in stats.keys():
        ax1.bar([title],stats["thickness"])
        ax2.bar([title],stats["thermal"])
        ax3.bar([title],stats["gprime"])
        ax4.bar([title],stats["ices"])
        ax1.tick_params(axis='x', labelrotation=45)
        ax2.tick_params(axis='x', labelrotation=45)
        ax3.tick_params(axis='x', labelrotation=45)
        ax4.tick_params(axis='x', labelrotation=45)

        ax1.set_title("H term")
        ax2.set_title("Thermal term")
        ax3.set_title("g' term")
        ax4.set_title("slope term")
