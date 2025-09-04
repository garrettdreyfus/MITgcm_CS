import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import xarray as xr
import gsw
import glob
from sklearn.linear_model import LinearRegression
from jmd95 import dens
from analysis import FStheory,slope,timeSeries, bottomMask, icemask, depthFromdZ 
from datainput import  matVarsFile, getIterNums, grabDeltaT, outPath, grabMatVars
from xmitgcm import open_mdsdataset
from matlabglib import GLIBfromFile
import matplotlib as mpl
from matplotlib.patches import Polygon
from matplotlib.animation import FFMpegFileWriter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from astropy.convolution import convolve, Box2DKernel
from tqdm import tqdm
from scipy.stats import pearsonr
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import cmocean
import os
from matplotlib.transforms import Affine2D
import ipdb



def barotropic_streamfunction_max(fname,times=np.array([]),res=1):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    shortname, fpath = outPath(fname) 
    depth = depthFromdZ(ds)
    UFULL = ds.UVEL.values
    VFULL = ds.VVEL.values
    bts = []
    for k in tqdm(range(0,ds.UVEL.values.shape[0],res)):
        U = UFULL[k,:,:,:]
        V = VFULL[k,:,:,:]
        fDZ = list(np.diff(ds.Z))
        fDZ.append(fDZ[-1])
        fDZ = np.asarray(fDZ)
        DZ = np.repeat(fDZ,U.shape[1]*U.shape[2])
        DZ = DZ.reshape(U.shape[0],U.shape[1],U.shape[2])
        UU = np.sum(U*DZ*ds.hFacW.values,axis=0)
        xs = ds.XC.values
        ys = list(np.diff(ds.YC.values))
        ys.append(ys[-1])
        ys = np.repeat(ys,U.shape[2])

        ys = ys.reshape(U.shape[1],U.shape[2],order="F")
        bt = np.cumsum(UU*ys,axis=0)
        bt[np.sum(ds.hFacC,axis=0)==0] = np.nan
        bts.append(np.nanmax(bt))
    return bts


def meltmapmovie(fname,description,times=np.array([])):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    shortname, fpath = outPath(fname) 
    moviewriter = FFMpegFileWriter(fps=1)
    fig,ax1 = plt.subplots()
    mask = np.logical_and(ds.hFacC.values[0]==0,np.sum(ds.hFacC.values,axis=0)!=0)
    theta = ds.THETA.values
    with moviewriter.saving(fig, fpath+"-meltmap.mp4" , dpi=100):
        for k in tqdm(range(int(ds.UVEL.values.shape[0]))):
            thetak = theta[k]
            thetaavg = np.nanmean(thetak,axis=0)
            melt= ds.SHIfwFlx.values[k]
            melt[~mask]=np.nan
            thetaavg[~mask]=np.nan
            frame = ax1.pcolormesh(-melt,cmap="jet",vmax=0.0012)
            if k==0:
                cb = plt.colorbar(frame)
            moviewriter.grab_frame()
            #cb.remove()
            frame.remove()

def bottomVtop(fname,description,times=np.array([])):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,prefix=["SALT","THETA","UVEL","VVEL","WVEL"],ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)


    inputname = "/".join(fname.split("/")[:-1])
    variables = grabMatVars(inputname,("h"))
    h = np.asarray(variables["h"]).T
    bmask = bottomMask(inputname,ds)
    icem = icemask(inputname,ds)
    shortname, fpath = outPath(fname) 
    fig,(ax1,ax2) = plt.subplots(1,2)
    theta = ds.THETA.values
    salt = ds.SALT.values
    uvel = ds.UVEL.values
    vvel = ds.WVEL.values
    X,Y = np.meshgrid(range(uvel.shape[3]),range(uvel.shape[2]))
    thetabot,saltbot,ubot,vbot = np.full_like(theta[1,1],0),np.full_like(theta[1,1],0),np.full_like(theta[1,1],0),np.full_like(theta[1,1],0)
    thetatop,salttop,utop,vtop = np.full_like(theta[1,1],0),np.full_like(theta[1,1],0),np.full_like(theta[1,1],0),np.full_like(theta[1,1],0)
    count=0
    for k in tqdm(range(int(ds.UVEL.values.shape[0]))):
        if times[k]>5:
            count+=1
            thetak = theta[k]
            thetabotc = thetak *bmask
            thetabotc[~bmask] = np.nan
            thetabot += np.nanmean(thetabotc,axis=0)

            saltk = salt[k]
            saltbotc = saltk *bmask
            saltbotc[~bmask] = np.nan
            #plt.imshow(saltbotc[:,:,90])
            #plt.show()
            #plt.imshow(ds.hFacC[:,:,90])
            #plt.show()
            saltbot += np.nanmean(saltbotc,axis=0)

            ubot += np.nanmean(uvel[k] * bmask,axis=0)
            vbot +=  np.nanmean(vvel[k] * bmask,axis=0)
            thetatopc = thetak *icem
            thetatopc[~icem] = np.nan
            thetatop += np.nanmean(thetatopc,axis=0)
            utop +=  np.nanmean(uvel[k] * icem,axis=0)
            vtop +=  np.nanmean(vvel[k] * icem,axis=0)
    thetabot = thetabot/count
    thetatop = thetatop/count
    saltbot = saltbot/count

    cax = ax1.pcolormesh(X,Y,thetatop,vmin=-1.9,vmax=-1.8,cmap=cmocean.cm.thermal)
    ax1.contour(X,Y,h,levels=range(-700,-500,25),colors="white",linestyles="solid")
    plt.colorbar(cax,ax=ax1)


    #cax = ax2.pcolormesh(X,Y,vbot,vmin=-0.0001,vmax=0.0001,cmap="jet")
    cax = ax2.pcolormesh(X,Y,thetabot,vmin=-2,vmax=-1.5,cmap=cmocean.cm.thermal)
    plt.colorbar(cax,ax=ax2)
    
    ax1.set_title("Top temp")
    ax2.set_title("Bottom temp")

    if k==1:
        cb = plt.colorbar(frame,ax=ax2,pad=0.05)

    plt.show()

def bottomspinup(fname,description,times=np.array([])):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,prefix=["SALT","THETA","UVEL","VVEL","WVEL","RHOAnoma"],ignore_unknown_vars=True,extra_variables = extra_variables)
    inputname = "/".join(fname.split("/")[:-1])
    bmask = bottomMask(inputname,ds)
    icem = icemask(inputname,ds)
    shortname, fpath = outPath(fname) 
    theta = ds.THETA.values
    salt = ds.SALT.values
    uvel = ds.UVEL.values
    vvel = ds.WVEL.values
    X,Y = np.meshgrid(range(uvel.shape[3]),range(uvel.shape[2]))
    thetabot,saltbot,ubot,vbot = np.full_like(theta[1,1],0),np.full_like(theta[1,1],0),np.full_like(theta[1,1],0),np.full_like(theta[1,1],0)
    thetatop,salttop,utop,vtop = np.full_like(theta[1,1],0),np.full_like(theta[1,1],0),np.full_like(theta[1,1],0),np.full_like(theta[1,1],0)
    count=0
    for k in tqdm(range(10)):
        fig,(ax1,ax2) = plt.subplots(1,2)
        thetak = theta[k]
        thetak[~bmask]=np.nan
        thetabot = np.nanmean(thetak,axis=0)
        saltk = salt[k]
        saltk[~bmask]=np.nan
        saltbot = np.nanmean(saltk,axis=0)
        uvelk = uvel[k]
        uvelk[ds.hFacC==0]=np.nan
        ubot = np.nansum(uvelk,axis=0)
        vvelk = vvel[k]
        vvelk[ds.hFacC==0]=np.nan
        vbot = np.nansum(np.abs(vvelk),axis=0)
        cax = ax1.pcolormesh(X,Y,thetabot,vmin=-1,vmax=1,cmap=cmocean.cm.thermal)
        plt.colorbar(cax,ax=ax1)
        cax = ax2.pcolormesh(X,Y,dens(saltbot,thetabot,400),cmap=cmocean.cm.haline,vmin=1029.50,vmax=1029.75)
        #cax = ax2.pcolormesh(X,Y,vbot,cmap=cmocean.cm.haline,vmin=-0.005,vmax=0.005)
        plt.colorbar(cax,ax=ax2)
        fig.canvas.manager.full_screen_toggle() # toggle fullscreen mode
        plt.show()


    ax1.set_title("Bottom temp")
    ax2.set_title("Bottom V")

    if k==1:
        cb = plt.colorbar(frame,ax=ax2,pad=0.05)

    plt.show()

def gprimeWidth(fname,xval,fig,ax1,title="",color="blue",marker="o"):
    data = timeSeries(fname)
    for k in data.keys():
        if k != "ts":
            try:
                data[k] = np.nanmean(data[k][data["ts"]>2])
            except:
                data[k]=np.nan
    if "gprime" not in data.keys():
        data["gprime"] = np.nan
    print(data["gprime"])
    variables = grabMatVars(fname,("Xeast","Xwest","Yicefront"))
    yice = float(variables["Yicefront"])
    shelf_width = float(variables["Xeast"])-float(variables["Xwest"])
    shiflx = -data["shiflx"]
    xval,shiflx = FStheory(fname,xval)
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))

    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    inputname = "/".join(fname.split("/")[:-1])
    icem = icemask(fname,ds)
    variables = grabMatVars(fname,("Hshelf","Xeast","Xwest","randtopog_height","Yicefront","Zcdw_pt_South","tEast","sEast","icedraft","zz","h","yy","xx"))
    icedraft = np.asarray(variables["icedraft"])
    h = np.asarray(variables["h"])
    tNorth = np.asarray(variables["tEast"])[-1,:]
    sNorth = np.asarray(variables["sEast"])[-1,:]
    rho0 = dens(sNorth,tNorth,0)
    zgl = np.nanmin(icedraft)

    max_height = variables["Zcdw_pt_South"][0][0]
    tcline_height = (max_height-75)/2.0+75
    localdens = dens(sNorth,tNorth,abs(200))
    d = localdens
    Zfull = np.asarray(list(ds.Z))
    if np.sum(abs(tNorth-tNorth[0])>0.2)>0:#and np.sum(t>0.5)>0:
        t = tNorth[sNorth>0.1]
        Z = Zfull[sNorth>0.1]
        s = sNorth[sNorth>0.1]
        mldi = np.where(abs(tNorth-tNorth[0])>0.2)[0][0]
        d = dens(sNorth,tNorth,Z[mldi])
        rho_1 = np.nanmean(d[:mldi])
        rho_2 = np.nanmean(d[mldi:min(mldi*2,len(d)-1)])
        gprime_ext = 9.8*(rho_2-rho_1)/np.mean((rho_1,rho_2))

    zz = np.asarray(variables["zz"])[0]
    zpyci = np.argmin(np.abs(np.abs(zz)-abs(max_height)))
    f = 1.3*10**-4
    ax1.scatter(gprime_ext*shelf_width,float(data["gprime"]),marker=marker,c=color)
    plt.xlabel(r'$g^{\prime}_{ext} (m/s^2)$',fontsize=18)
    plt.ylabel(r'$g^{\prime}_{in} (m/s^2)$',fontsize=18)
    return gprime_ext

def connectionPlot(fname,xval,fig,ax1,title="",color="blue",marker="o"):
    data = timeSeries(fname)
    variables = grabMatVars(fname,("Xeast","Xwest","Yicefront","saltfluxvals"))
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    dx = np.gradient(ds.XG)[0]
    dy = np.gradient(ds.YG)[0] 
    saltfluxvals = np.asarray(variables["saltfluxvals"])
    meltsaltflux = np.sum(ds.SHIfwFlx,axis=[1,2]).values*dx*dy*34.5
    polynaflux = np.sum(saltfluxvals*dx*dy)
    starttime=7
    ax1.scatter(np.mean(-meltsaltflux[data["ts"]>starttime]+polynaflux),np.mean(data["overturning_connect"][data["ts"]>starttime]),marker=marker,c=color,label=title,s=100)

    return 0

def gprimeAll(fname,xval,fig,ax1,title="",color="blue",marker="o"):
    data = timeSeries(fname)
    for k in data.keys():
        if k != "ts":
            try:
                data[k] = np.nanmean(data[k][data["ts"]>2])
            except:
                data[k]=np.nan
    if "gprime" not in data.keys():
        data["gprime"] = np.nan
    variables = grabMatVars(fname,("Xeast","Xwest","Yicefront"))
    yice = float(variables["Yicefront"])
    shelf_width = float(variables["Xeast"])-float(variables["Xwest"])
    shiflx = -data["shiflx"]
    xval,shiflx = FStheory(fname,xval)
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    extra_variables["KPPdiffT"] = dict(dims=["k","j","i"], attrs=dict(standard_name="KPP mld", units="m"))

    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    inputname = "/".join(fname.split("/")[:-1])
    icem = icemask(fname,ds)
    variables = grabMatVars(fname,("Hshelf","Xeast","Xwest","randtopog_height","Yicefront","Zcdw_pt_South","tEast","sEast","icedraft","zz","h","yy","xx"))
    icedraft = np.asarray(variables["icedraft"])
    h = np.asarray(variables["h"])
    #
    tNorth = np.asarray(variables["tEast"])[-1,:]
    sNorth = np.asarray(variables["sEast"])[-1,:]
    rho0 = dens(sNorth,tNorth,0)
    zgl = np.nanmin(icedraft)

    max_height = variables["Zcdw_pt_South"][0][0]
    tcline_height = (max_height-75)/2.0+75
    localdens = dens(sNorth,tNorth,abs(tcline_height))
    d = localdens
    Z = np.asarray(list(ds.Z))
    zz = np.asarray(variables["zz"])[0]
    if np.sum(abs(tNorth-tNorth[0])>0.2)>0:#and np.sum(t>0.5)>0:
        mldi = np.where(abs(tNorth-tNorth[0])>0.2)[0][0]
        #cdwi = np.where(t>0)[0][0]
        rho_1 = np.nanmean(d[:mldi])
        rho_2 = np.nanmean(d[mldi:min(mldi*2,len(d)-1)])
        #gprime_ext = np.nanmean(np.abs(np.diff(d)/np.diff(zz)))
        gprime_ext = 9.8*(rho_2-rho_1)/np.mean((rho_1,rho_2))

    zz = np.asarray(variables["zz"])[0]
    zpyci = np.argmin(np.abs(np.abs(zz)-abs(max_height)))
    gprime_ext = 9.8*(np.mean(localdens[zpyci:min(zpyci*2,len(localdens)-1)])-np.mean(localdens[:zpyci]))/np.mean(localdens[:min(zpyci*2,len(localdens)-1)])
    f = 1.3*10**-4
    #ax1.scatter(gprime_ext*shelf_width,float(data["gprime"])*840,marker=marker,c=color)
    #ax1.scatter(gprime_ext*shelf_width,float(data["gprime"])*625000-525,marker=marker,c=color)
    #ax1.scatter(gprime_ext*shelf_width,float(data["gprime"]),marker=marker,c=color)
    W0 = gprime_ext*shelf_width/float(data["gprime"])
    kpp = np.nanmean(ds.KPPdiffT.values,axis=0)
    kpp = np.nanmean(kpp,axis=0)
    kpp = np.nanmean(kpp[ds.YC.values<60000])
    small=np.log10(np.abs(ds.THETA_inst.values[1]))<-3
    smallmask=np.logical_and(small,ds.hFacC!=0)
    smallmask = np.sum(smallmask,axis=0)
    ax1.scatter(float(data["gprime"]),np.sum(smallmask[YY<60000]),marker=marker,c=color)
    plt.xlabel(r'$g^{\prime}_{ext} (m/s^2)$',fontsize=18)
    plt.ylabel(r'$g^{\prime}_{in} (m/s^2)$',fontsize=18)
    return gprime_ext

def saltBudget(fname,xval,fig,ax1,title="",color="blue",marker="o"):
    data = timeSeries(fname)
    if "tcdw" not in data:
        return 0
    for k in data.keys():
        if k != "ts":
            try:
                data[k] = np.nanmean(data[k][data["ts"]>2])
            except:
                data[k]=np.nan
    if "gprime" not in data.keys():
        data["gprime"] = np.nan
    variables = grabMatVars(fname,("Xeast","Xwest","Yicefront"))
    yice = float(variables["Yicefront"])
    shelf_width = float(variables["Xeast"])-float(variables["Xwest"])
    shiflx = -data["shiflx"]
    xval,shiflx = FStheory(fname,xval)
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))

    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    minvel = (ds.UVEL*(ds.THETA+1.8)).mean(dim="time",skipna=True).min(dim="XG").min(dim="XC").min(dim="Z")#.argmin(dim="YC")
    ht = ds.UVEL.values*(ds.THETA.values+1.8)
    ht = np.nanmean(ht,axis=0)
    ht = np.nanmin(ht,axis=0)
    ht = np.mean(ht,axis=1)
    minvel=ht

    inputname = "/".join(fname.split("/")[:-1])
    icem = icemask(fname,ds)
    variables = grabMatVars(fname,("Hshelf","Xeast","Xwest","randtopog_height","Yicefront","Zcdw_pt_South","tEast","sEast","icedraft","zz","h","yy","xx"))
    icedraft = np.asarray(variables["icedraft"])
    h = np.asarray(variables["h"])
    meltline = ds.SHIfwFlx.mean(dim="time").mean(dim="XC").values
    meltline = minvel
    weighted = np.sum(meltline*np.asarray(range(len(meltline))))/np.sum(range(len(meltline)))
    #iceline = np.mean(icedraft,axis=0)
    #print("velmin: ",ds.UVEL.mean(dim="time",skipna=True).min(dim="XG").min(dim="Z").argmin(dim="YC").values)
    #plt.plot(iceline)
    #ax2 = plt.gca().twinx()
    #plt.title(fname)
    #ax2.plot(meltline)
    #plt.show()
    #print(minvel)
    #
    tNorth = np.asarray(variables["tEast"])[-1,:]
    sNorth = np.asarray(variables["sEast"])[-1,:]
    rho0 = dens(sNorth,tNorth,0)
    zgl = np.nanmin(icedraft)

    max_height = variables["Zcdw_pt_South"][0][0]
    tcline_height = (max_height-75)/2.0+75
    localdens = dens(sNorth,tNorth,abs(tcline_height))
    d = localdens
    Z = np.asarray(list(ds.Z))
    zz = np.asarray(variables["zz"])[0]
    if np.sum(abs(tNorth-tNorth[0])>0.2)>0:#and np.sum(t>0.5)>0:
        mldi = np.where(abs(tNorth-tNorth[0])>0.2)[0][0]
        #cdwi = np.where(t>0)[0][0]
        rho_1 = np.nanmean(d[:mldi])
        rho_2 = np.nanmean(d[mldi:min(mldi*2,len(d)-1)])
        #gprime_ext = np.nanmean(np.abs(np.diff(d)/np.diff(zz)))
        gprime_ext = 9.8*(rho_2-rho_1)/np.mean((rho_1,rho_2))

    zz = np.asarray(variables["zz"])[0]
    zpyci = np.argmin(np.abs(np.abs(zz)-abs(max_height)))
    gprime_ext = 9.8*(np.mean(localdens[zpyci:min(zpyci*2,len(localdens)-1)])-np.mean(localdens[:zpyci]))/np.mean(localdens[:min(zpyci*2,len(localdens)-1)])
    f = 1.3*10**-4
    shelf_width = float(variables["Xeast"])-float(variables["Xwest"])
    #ax1.scatter(float(data["scdw"]-data["ssurf"])/float(data["tcdw"]-data["tsurf"]),data["froude"],marker=marker,c=color)
    ax1.scatter(shelf_width,float(data["tcdw"]-data["tsurf"])**2,marker=marker,c=color)
    #print(sNorth)
    #ax1.scatter(float(data["scdw"]-data["ssurf"])/float(data["tcdw"]-data["tsurf"]),np.std(sNorth),marker=marker,c=color)
    plt.xlabel(r'Width (m)',fontsize=18)
    plt.ylabel(r'$(T_{cdw} - T_{surf})^2$',fontsize=18)
    return gprime_ext

def steadyStateAverageSimple(fname,xval,fig,ax1,title="",color="blue",marker="o"):
    shortname,_ = outPath(fname)
    shortname = shortname.split("|")[1]
    if xval == None:
        xval,shiflx,stats = FStheory(fname,xval)
        if ~np.isnan(np.asarray(xval)).all():
                ax1.scatter(xval,shiflx,c=color,marker=marker,label=shortname,s=125)
                # ax1.plot(xval,shiflx,c=color,marker=marker,label=shortname)
        return xval

    data = timeSeries(fname)
    print("XVAL SUPPLIED")
    print(data["ts"])
    shiflx = -np.nanmean(data["shiflx"][data["ts"]>7])
    ax1.scatter(xval,shiflx,c=color,marker=marker,label=shortname,s=125)
    return xval

def steadyStateHT(fname,xval,fig,ax1,title="",color="blue",marker="o"):
    data = timeSeries(fname)
    for k in data.keys():
        if k != "ts":
            try:
                data[k] = np.nanmean(data[k][data["ts"]>2])
            except:
                data[k]=np.nan
    if "gprime" not in data.keys():
        data["gprime"] = np.nan
    variables = grabMatVars(fname,("Xeast","Xwest","Yicefront"))
    yice = float(variables["Yicefront"])
    shelf_width = float(variables["Xeast"])-float(variables["Xwest"])
    shiflx = -data["shiflx"]
    xval,shiflx = FStheory(fname,xval)
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))

    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    inputname = "/".join(fname.split("/")[:-1])
    icem = icemask(fname,ds)

    ax1.scatter(xval,-float(data["incavity"])/np.sum(icem),c=color,marker=marker)
    #ax1.scatter(xval,np.sum(icem),c=color,marker=marker)
    return data["incavity"]


def timeSeriesDashboard(fname,label,fig,axises,times=np.array([]),color="yellow"):
    ((ax1,ax2,ax5),(ax3,ax4,ax6)) = axises 
    data = timeSeries(fname)
    starttime = 0

    ds = open_mdsdataset(fname,prefix=["SHIfwFlx"])

    dx = np.gradient(ds.XG)[0]
    dy = np.gradient(ds.YG)[0] 

    variables = grabMatVars(fname,("saltfluxvals"))
    saltfluxvals = np.asarray(variables["saltfluxvals"])

    meltsaltflux = np.sum(ds.SHIfwFlx,axis=[1,2]).values*dx*dy*34.5
    polynaflux = np.sum(saltfluxvals*dx*dy)

    variables = grabMatVars(fname,("tEast",'sEast','zz'))
    tNorth = np.asarray(variables["tEast"])[-1,:]
    sNorth = np.asarray(variables["sEast"])[-1,:]
    localdens = dens(sNorth,tNorth,abs(200))
    zz = np.asarray(variables["zz"])[0]
    d = localdens
    Zfull = np.asarray(list(ds.Z))
    gradd = np.abs(np.diff(localdens)/np.diff(zz))
    #average depth of all above 80th percentile
    tcline_height=np.mean(zz[:-1][gradd>np.quantile(gradd,0.85)])#+75
    zpyci = np.argmin(np.abs(np.abs(zz)-abs(tcline_height)))
    localdens = dens(sNorth,tNorth,abs(zz[zpyci]))

    ## calculation of gprime
    rho_1i = np.logical_and(zz<zz[zpyci],zz>zz[zpyci]-50)
    rho_2i = np.logical_and(zz<zz[zpyci]+50,zz>zz[zpyci])
    gprime_ext = 9.8*(np.nanmean(localdens[rho_1i])-np.nanmean(localdens[rho_2i]))/np.mean(localdens[np.logical_or(rho_1i,rho_2i)])


    if len(data["ts"])==0:
        print(fname, "EMPTY")
        return 
    if np.nanmax(data["ts"][~np.isnan(data["theta"])])<9:
        print(fname,np.nanmax(data["ts"][~np.isnan(data["theta"])]))


    ax1.plot(data["ts"][data["ts"]>starttime],data["theta"][data["ts"]>starttime],label=label,c=color)
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Potential Temperature")
    ## salt plot
    ax2.plot(data["ts"][data["ts"]>starttime],data["salt"][data["ts"]>starttime],label=label,c=color)
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Salinity")

    ## kinetic energy plot
    kes = []
    #ax3.scatter(gprime_ext*np.mean(meltsaltflux[1:]/polynaflux),np.mean(data["offshorefraction"][data["ts"]>starttime][6:]),label=label,c=color)
    #ax3.scatter(np.mean(meltsaltflux[1:]/polynaflux),np.mean(data["offshorefraction"][data["ts"]>starttime][6:]),label=label,c=color)
    ax3.scatter(np.mean(-meltsaltflux[1:][data["ts"]>starttime]+polynaflux),np.mean(data["overturning_connect"][data["ts"]>starttime][starttime:]),label=label,c=color)
    ax3.set_xlabel("$B_{total}$")
    ax3.set_ylabel("Fractional streamline connectedness")

    ax4.plot(data["ts"][data["ts"]>starttime],data["shiflx"][data["ts"]>starttime],c=color)

    ax4.set_xlabel("Time")
    ax4.set_ylabel("Melt Rate m/yr")

    #ax5.plot(data["ts"][data["ts"]>starttime],data["polynasaltdiff"][data["ts"]>starttime],c=color)
    ax5.plot(data["overturning_connect"][data["ts"]>starttime],meltsaltflux[1:][data["ts"]>starttime]-polynaflux,c=color)
    ax5.set_xlabel("overturning")
    ax5.set_ylabel("Btotal")

    #ax6.plot(data["ts"][data["ts"]>starttime],data["incavity"][data["ts"]>starttime],c=color)
    ax6.plot(data["ts"][data["ts"]>starttime],data["overturning_connect"][data["ts"]>starttime],c=color)
    ax6.set_xlabel("Time")
    ax6.set_ylabel("Meltflux/polynaflux")

    if np.max(data["ts"][~np.isnan(data["shiflx"])])<7.5:
        print(label)
    
    ## kinetic energy plot
    #totalSHIfwFlx = (ds.SHIfwFlx*ds.XC*ds.YC*ds.Z*ds.hFacC).sum(axis=[1,2,3]).values
    bottomtemps = []
    surfacetemps = []
 
def crossSectionAverage(fname,description,selval,quant="THETA",dim="zonal",ax1=None,show=True,savepath=False,fixcb=True):
    if not ax1:
        fig,ax1 = plt.subplots(figsize=(10,8))
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    if quant=="UVEL":
        ds = open_mdsdataset(fname,prefix=["UVEL"],ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
        ds = ds.mean(dim="time")
        if dim == "meridional":
            ind = np.argmin(np.abs(ds.XC.values-selval))
            zonal_average = ds.isel(XG=ind,drop=True)
            zonal_average = zonal_average.isel(XC=ind,drop=True)
            ys = zonal_average.YC.values
        if dim == "zonal":
            ind = np.argmin(np.abs(ds.YC.values-selval))
            zonal_average = ds.isel(YG=ind,drop=True)
            zonal_average = zonal_average.isel(YC=ind,drop=True)
            ys = zonal_average.XC.values
        zvals = zonal_average[quant].values
        zvals[zonal_average.hFacC.values==0] = np.nan
        m = np.nanmedian(zvals)
        s = np.nanstd(zvals)
        tmin, tmax = m-2*s,m+s*2
        zs = zonal_average.Z.values
        length = zvals.shape[0]
        shortname, fpath = outPath(fname) 
        fig.suptitle(description)
    elif quant!="DENS":
        ds = open_mdsdataset(fname,prefix=["THETA","SALT","WVEL","VVEL","UVEL","RHOAnoma"],ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
        ds = ds.mean(dim="time")
        if dim == "meridional":
            ind = np.argmin(np.abs(ds.XC.values-selval))
            zonal_average = ds.isel(XC=ind,drop=True)
            ys = zonal_average.YC.values
        if dim == "zonal":
            ind = np.argmin(np.abs(ds.YC.values-selval))
            zonal_average = ds.isel(YC=ind,drop=True)
            ys = zonal_average.XC.values
        zvals = zonal_average[quant].values
        zvals[zonal_average.hFacC.values==0] = np.nan
        m = np.nanmedian(zvals)
        s = np.nanstd(zvals)
        tmin, tmax = m-2*s,m+s*2
        zs = zonal_average.Z.values
        length = zvals.shape[0]
        shortname, fpath = outPath(fname) 
        fig.suptitle(description)
    if quant=="DENS":
        ds = open_mdsdataset(fname,prefix=["SALT","THETA"],ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
        ds = ds.mean(dim="time")
        #zonal_average = ds.where(ds.hFacC==1).isel(XC=100)
        if dim == "meridional":
            ind = np.argmin(np.abs(ds.XC.values-selval))
            #zonal_average = ds.where(ds.hFacC==1).mean(dim="XC",skipna=True)
            zonal_average = ds.isel(XC=ind)
            ys = zonal_average.YC.values
        if dim == "zonal":
            ind = np.argmin(np.abs(ds.YC.values-selval))
            zonal_average = ds.isel(YC=ind)
            ys = zonal_average.XC.values
        shortname, fpath = outPath(fname) 
        fig.suptitle(description)
        zs = zonal_average.Z.values
        #tmin, tmax = np.nanmin(zonal_average[quant]), np.nanmax(zonal_average[quant])
        zvals = (zonal_average["SALT"].values,zonal_average["THETA"].values)
        length = zvals[0].shape[0]
        zvals = gsw.sigma0(zvals[0],zvals[1])
        zvals[zonal_average.hFacC.values==0] = np.nan

    if fixcb:
        if quant=="THETA":
            c=ax1.pcolormesh(ys/1000,zs/1000,zvals,cmap=cmocean.cm.thermal,zorder=5,vmin=-2,vmax=-1.8)
        if quant=="SALT":
            c=ax1.pcolormesh(ys/1000,zs/1000,zvals,cmap=cmocean.cm.thermal,zorder=5,vmin=34.1,vmax=34.4)
        if "VEL" in quant:
            c=ax1.pcolormesh(ys/1000,zs/1000,zvals,cmap="seismic",zorder=5,vmin=-0.075,vmax=0.075)
        if quant=="DENS":
            c = ax1.pcolormesh(ys,zs,zvals,cmap="jet",vmin=27.4,vmax=27.7)
    else:
        variables = grabMatVars(fname,("Yicefront"))
        yice = float(variables["Yicefront"])
        cavmin = np.nanmin(zvals[:,ys<yice])
        cavmax = np.nanmax(zvals[:,ys<yice])
        if quant=="THETA":
            c=ax1.pcolormesh(ys/1000,zs/1000,zvals,cmap=cmocean.cm.thermal,zorder=5,vmin=cavmin,vmax=cavmax)
        if quant=="SALT":
            c=ax1.pcolormesh(ys/1000,zs/1000,zvals,cmap=cmocean.cm.thermal,zorder=5,vmin=cavmin,vmax=cavmax)
        if "VEL" in quant:
            c=ax1.pcolormesh(ys/1000,zs/1000,zvals,cmap="seismic",zorder=5,vmin=cavmin,vmax=cavmax)
        if quant=="DENS":
            c = ax1.pcolormesh(ys,zs,zvals,cmap="jet",vmin=cavmin,vmax=cavmax)
        else:
            c = ax1.pcolormesh(ys,zs,zvals,cmap="jet",vmin=cavmin,vmax=cavmax)
 
        
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    variables = grabMatVars(fname,("Xeast","Xwest","Yicefront"))
    yice = float(variables["Yicefront"])/1000
    ices = slope(fname)
    print("slope: ",ices)

        
    caxout = inset_axes(
        ax1,
        width="2%",  # width: 5% of parent_bbox width
        height="40%",  # height: 50%
        loc="lower left",
        bbox_to_anchor=(1.0, 0.30, 1, 1),
        bbox_transform=ax1.transAxes,
        borderpad=1,
    )
    caxout.tick_params(labelsize=18)
    plt.colorbar(c,cax=caxout)#,ticks=[-2,0,1])
    caxout.set_ylabel('$^\circ$C', rotation=0,fontsize=18)

    ax1.set_ylabel('Depth (km)',fontsize=18)
    ax1.set_xlabel('Cross Shelf Distance (km)',fontsize=18)
    if show:
        plt.show()
    if savepath:
        plt.savefig(savepath)
    plt.close()


def simpleSurfaceAverage(fname,description,quant="THETA",res=1,dim="zonal",ax1=None,show=True):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,prefix=["THETA","SALT","WVEL","VVEL"],ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    ds = ds.mean(dim="time")
    #plt.plot(range(len(ds.Z.values)),ds.Z.values)
    #plt.show()
    depth_slice = ds.isel(Z=20,drop=True)
    theta = depth_slice.THETA.values
    theta[depth_slice.hFacC.values==0] = np.nan
    salt = depth_slice.SALT.values
    salt[depth_slice.hFacC.values==0] = np.nan
    fig, (ax1,ax2) = plt.subplots(1,2)
    im = ax1.imshow(theta,vmin=-2,vmax=-1.8,cmap=cmocean.cm.thermal)
    plt.colorbar(im,ax=ax1,orientation="horizontal")
    im = ax2.imshow(salt,vmin=34.2,vmax=34.4,cmap=cmocean.cm.haline)
    plt.colorbar(im,ax=ax2,orientation="horizontal")
    plt.show()



def barotropic_streamfunction_graph(fname,description,times=np.array([]),res=1):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    #times=getIterNums(fname)
    #ds = open_mdsdataset(fname,prefix=["UVEL","VVEL"],ignore_unknown_vars=True)#,iters=times)
    ds = open_mdsdataset(fname,prefix=["SALT","THETA","UVEL","VVEL","WVEL"],ignore_unknown_vars=True,extra_variables = extra_variables)
    shortname, fpath = outPath(fname) 
    depth = depthFromdZ(ds)
    UFULL = ds.UVEL.values#
    VFULL = ds.VVEL.values
    #cavity_mask = (ds.SHIfwFlx[0].values !=0)
    #for k in tqdm(range(int(ds.UVEL.values.shape[0]*3/4),ds.UVEL.values.shape[0],res)):
    for k in tqdm(range(10)):
        fig,ax1 = plt.subplots()
        U = UFULL[k]
        V = VFULL[k]
        fDZ = list(np.diff(ds.Z))
        fDZ.append(fDZ[-1])
        fDZ = np.asarray(fDZ)
        DZ = np.repeat(fDZ,U.shape[1]*U.shape[2])
        DZ = DZ.reshape(U.shape[0],U.shape[1],U.shape[2])
        UU = np.sum(U*DZ*ds.hFacW.values,axis=0)
        xs = ds.XC.values
        ys = list(np.diff(ds.YC.values))
        ys.append(ys[-1])
        ys = np.repeat(ys,U.shape[2])
        ys = ys.reshape(U.shape[1],U.shape[2],order="F")
        bt = np.cumsum(UU*ys,axis=0)
        bt[np.sum(ds.hFacC,axis=0)==0] = np.nan
        m,s = np.nanmedian(bt[ds.YC.values<60000]),np.nanstd(bt[ds.YC.values<60000])
        #vmin,vmax = m-3*s,m+3*s
        vmin,vmax=0,30000
        frame = ax1.pcolormesh(ds.XC.values,ds.YC.values,bt,vmin=vmin,vmax=vmax)
        ax1.contour(ds.XC.values,ds.YC.values,bt,levels=60,colors="black",vmin=vmin,vmax=vmax)
        ax1.quiver(ds.XC.values,ds.YC.values,np.sum(U,axis=0),np.sum(V,axis=0))
        plt.title(str(k))
        cb = plt.colorbar(frame)
        plt.show()



def barotropic_streamfunction_max(fname,times=np.array([]),res=1):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    shortname, fpath = outPath(fname) 
    depth = depthFromdZ(ds)
    UFULL = ds.UVEL.values
    VFULL = ds.VVEL.values
    bts = []
    for k in tqdm(range(0,ds.UVEL.values.shape[0],res)):
        U = UFULL[k,:,:,:]
        V = VFULL[k,:,:,:]
        fDZ = list(np.diff(ds.Z))
        fDZ.append(fDZ[-1])
        fDZ = np.asarray(fDZ)
        DZ = np.repeat(fDZ,U.shape[1]*U.shape[2])
        DZ = DZ.reshape(U.shape[0],U.shape[1],U.shape[2])
        UU = np.sum(U*DZ*ds.hFacW.values,axis=0)
        xs = ds.XC.values
        ys = list(np.diff(ds.YC.values))
        ys.append(ys[-1])
        ys = np.repeat(ys,U.shape[2])

        ys = ys.reshape(U.shape[1],U.shape[2],order="F")
        bt = np.cumsum(UU*ys,axis=0)
        bt[np.sum(ds.hFacC,axis=0)==0] = np.nan
        bts.append(np.nanmax(bt))
    return bts

def circulationFigure(fname,description,times=np.array([])):
    mpl.rcParams['savefig.dpi'] = 500
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,prefix=["THETA","SALT","momKE","SHIfwFlx","VVEL","UVEL","WVEL"],ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    fig,ax1 = plt.subplots()
    mask = np.logical_and(ds.hFacC.values[0]==0,np.sum(ds.hFacC.values,axis=0)!=0)
    mask[1,:]=0
    zs = ds.Z.values
    xs = ds.XG.values
    ys = ds.YG.values
    times = np.asarray(times)
    ts = np.asarray(times/60.0/60.0/24.0/365.0)
    ts = ds.time.values*grabDeltaT(fname)/60.0/60.0/24.0/365.0
    tsnew = np.full_like(ts,0,dtype=float)
    tsnew[:] = ts
    ts = tsnew
    ts = ts/1000000000

    glib = GLIBfromFile(matVarsFile(fname))

    uvel = np.mean(ds.UVEL.values[ts>5],axis=0)
    vvel = np.mean(ds.VVEL.values[ts>5],axis=0)

    theta= np.mean(ds.SALT.values[ts>5],axis=0)
    hfac = ds.hFacC.values
    interfacez = np.zeros(theta.shape[1:])

    interfaceu = np.zeros(theta.shape[1:])
    interfacev = np.zeros(theta.shape[1:])
    bottomz = np.zeros(theta.shape[1:])
    for i in range(theta.shape[1]):
        for j in range(theta.shape[2]):
            zc = np.where(np.diff(np.sign(theta[:,i,j]-34.3)))[0]
            if len(zc)>0:
                interfacez[i,j] = zs[zc[0]]
                interfaceu[i,j] = np.nanmean((uvel[zc[0]:,i,j]*zs[zc[0]:]))/np.sum(zs[zc[0]:])
                interfacev[i,j] = np.nanmean((vvel[zc[0]:,i,j]*zs[zc[0]:]))/np.sum(zs[zc[0]:])
            else:
                interfacez[i,j] = np.nan
                interfaceu[i,j] = np.nan
                interfacev[i,j] = np.nan
            if np.sum(hfac[:,i,j]!=0)>0:
                bottomz[i,j] = np.nanmin(zs[hfac[:,i,j]!=0])
            else:
                bottomz[i,j] = np.nan

            
    ax1.set_facecolor("#BBAF98")
    xs,ys=xs/1000,ys/1000
    X,Y= np.meshgrid(xs,ys)

    interfacez[np.logical_and(X>275,Y<200)]=np.nan
    bottomz[np.logical_and(X>275,Y<200)]=np.nan
    #plt.imshow(np.logical_and(X>275,Y<200))
    #plt.show()
 
    c= plt.pcolormesh(xs,ys,interfacez[:,::-1],cmap=cmocean.cm.deep)

    caxout = plt.colorbar(c,aspect=40,shrink=0.8,ticks=range(-900,-199,175))
    caxout.ax.tick_params(labelsize=18)
    #plt.quiver(X[::5,::5],Y[::5,::5],interfaceu[::5,::5],interfacev[::5,::5],color="white")
    #plt.quiver(X[::1,::1],Y[::1,::1],interfaceu[::1,::1],interfacev[::1,::1],color="white")
    box_kernel = Box2DKernel(5,mode="center")
    X = convolve(X, box_kernel,preserve_nan=True)
    Y = convolve(Y, box_kernel,preserve_nan=True)
    interfaceu = convolve(interfaceu, box_kernel,preserve_nan=True,boundary=None)
    interfacev = convolve(interfacev, box_kernel,preserve_nan=True,boundary=None)
    interfaceu[np.isnan(interfacez)]=np.nan
    interfaceu[interfaceu==0]=np.nan
    interfacev[np.isnan(interfacez)]=np.nan
    interfacev[interfacev==0]=np.nan

    plt.quiver(X[::5,::5],Y[::5,::5],interfaceu[::5,::5],interfacev[::5,::5],color="white",scale=0.0025,zorder=5,width=0.005)
    plt.hlines(150,126,271,linewidth=7.5,color="orange")

    plt.gca().tick_params(labelsize=15)
    glibi = np.argmin(np.abs(np.abs(zs)-np.abs(glib)))
    plt.contour(xs,ys,bottomz[:,::-1],[zs[glibi+1]],colors=["red"],linewidths=3)
    plt.ylabel(r'y (km)',fontsize=18)
    plt.xlabel(r'x (km)',fontsize=18)
    plt.tight_layout()
    plt.show()



def meltmapmovie(fname,description,times=np.array([])):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    shortname, fpath = outPath(fname) 
    moviewriter = FFMpegFileWriter(fps=1)
    fig,ax1 = plt.subplots()
    mask = np.logical_and(ds.hFacC.values[0]==0,np.sum(ds.hFacC.values,axis=0)!=0)
    theta = ds.THETA.values
    with moviewriter.saving(fig, fpath+"-meltmap.mp4" , dpi=100):
        for k in tqdm(range(int(ds.UVEL.values.shape[0]))):
            thetak = theta[k]
            thetaavg = np.nanmean(thetak,axis=0)
            melt= ds.SHIfwFlx.values[k]
            melt[~mask]=np.nan
            thetaavg[~mask]=np.nan
            frame = ax1.pcolormesh(-melt,cmap="jet",vmax=0.0012)
            if k==0:
                cb = plt.colorbar(frame)
            moviewriter.grab_frame()
            #cb.remove()
            frame.remove()

def mixmap(fname,description,times=np.array([])):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    shortname, fpath = outPath(fname) 
    moviewriter = FFMpegFileWriter(fps=1)
    fig,ax1 = plt.subplots()

    mask = np.logical_and(ds.hFacC.values[0]==0,np.sum(ds.hFacC.values,axis=0)!=0)
    theta = ds.THETA.values
    with moviewriter.saving(fig, fpath+"-mixmap.mp4" , dpi=100):
        for k in tqdm(range(int(ds.UVEL.values.shape[0]))):
            thetak = theta[k]
            thetaavg = np.nanmean(thetak,axis=0)
            melt= ds.MXLDEPTH.values[k]
            melt[~mask]=np.nan
            thetaavg[~mask]=np.nan
            frame = ax1.pcolormesh(melt,cmap="jet")
            if k==0:
                cb = plt.colorbar(frame)
            moviewriter.grab_frame()
            #cb.remove()
            frame.remove()


def topMap(fname,description,quant="THETA",times=np.array([]),show=False,savepath=False,fixcb=True):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))

    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    inputname = "/".join(fname.split("/")[:-1])
    bmask = bottomMask(inputname,ds)
    icem = icemask(inputname,ds)
    shortname, fpath = outPath(fname) 
    #fig,ax1 = plt.subplots(1,1)
    theta = ds.THETA.values
    X = ds.XC.values
    Y = ds.YC.values
    thetatop = np.full_like(theta[1,1],0)
    count=0
    #for k in tqdm(range(int(theta.shape[0]))):
        #if times[k]>5:
            #count+=1
            #thetak = theta[k]
            #thetatopc = thetak *icem
            #thetatopc[~icem]=np.nan
            ##thetatopc=thetak[0]
            #if np.sum(~np.isnan(thetatopc)) !=0:
                #thetatop += np.nanmean(thetatopc,axis=0)
    #thetatop = thetatop/count
    #mask = np.logical_and(ds.hFacC.values[0]==0,np.sum(ds.hFacC.values,axis=0)!=0)
    #thetatop[~mask]=np.nan
    #
    #cax = ax1.pcolormesh(X/1000,Y/1000,thetatop,cmap=cmocean.cm.thermal,vmin=-0.8,vmax=0.2)
    #cb = plt.colorbar(cax,ax=ax1,pad=0.05)
#
    #plt.xlabel("x (km)",fontsize=18)
    #plt.ylabel("y (km)",fontsize=18)
    #plt.xlim(50,350)
    #plt.ylim(0,150)
    #plt.xticks(fontsize=18)
    #plt.yticks(fontsize=18)
#
    ##ax1.set_title("Width: 100",fontsize=18)
    #plt.tight_layout
    #plt.show()


    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))

    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    inputname = "/".join(fname.split("/")[:-1])
    bmask = bottomMask(inputname,ds)
    icem = icemask(inputname,ds)
    shortname, fpath = outPath(fname) 
    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(18,12))
    theta = ds[quant].values
    uvel = ds.UVEL.values
    vvel = ds.VVEL.values
    X,Y = np.meshgrid(range(uvel.shape[3]),range(uvel.shape[2]))
    thetabot,ubot,vbot = np.full_like(theta[1,1],0),np.full_like(theta[1,1],0),np.full_like(theta[1,1],0)
    thetatop,utop,vtop = np.full_like(theta[1,1],0),np.full_like(theta[1,1],0),np.full_like(theta[1,1],0)
    count=0
    for k in tqdm(range(int(ds.UVEL.values.shape[0]))):
        if times[k]>5:
            count+=1
            thetak = theta[k]
            thetabotc = thetak *bmask
            thetabotc[~bmask] = np.nan
            thetabot += np.nanmean(thetabotc,axis=0)
            ubot += np.nanmean(uvel[k] * bmask,axis=0)
            vbot +=  np.nanmean(vvel[k] * bmask,axis=0)
            thetatopc = thetak *icem
            thetatopc[~icem] = np.nan
            thetatop += np.nanmean(thetatopc,axis=0)
            utop +=  np.nanmean(uvel[k] * icem,axis=0)
            vtop +=  np.nanmean(vvel[k] * icem,axis=0)
    thetabot = thetabot/count
    thetatop = thetatop/count

    variables = grabMatVars(fname[:-1],("Yicefront"))
    yice = float(variables["Yicefront"])
    if quant=="THETA":
        if fixcb:

                ax1.pcolormesh(X,Y,thetabot,vmin=-2.2,vmax=-1.8,cmap=cmocean.cm.thermal)
                frame = ax2.pcolormesh(X,Y,thetatop,vmin=-2.2,vmax=-1.8,cmap=cmocean.cm.thermal)
                cb = plt.colorbar(frame,ax=ax2,pad=0.05)
        else:
                cavmin = np.nanmin(thetabot[Y<yice])
                cavmax = np.nanmax(thetabot[Y<yice])
                frame = ax1.pcolormesh(X,Y,thetabot,cmap=cmocean.cm.thermal,vmin=cavmin,vmax=cavmax)
                cb = plt.colorbar(frame,ax=ax1,pad=0.05)
                cavmin = np.nanmin(thetatop[Y<yice])
                cavmax = np.nanmax(thetatop[Y<yice])
                frame = ax2.pcolormesh(X,Y,thetatop,cmap=cmocean.cm.thermal,vmin=cavmin,vmax=cavmax)
                cb = plt.colorbar(frame,ax=ax2,pad=0.05)
    if quant=="SALT":
        if fixcb:
                ax1.pcolormesh(X,Y,thetabot,vmin=34.15,vmax=34.5,cmap=cmocean.cm.haline)
                frame = ax2.pcolormesh(X,Y,thetatop,vmin=34.15,vmax=34.4,cmap=cmocean.cm.haline)
                cb = plt.colorbar(frame,ax=ax2,pad=0.05)
        else:
                cavmin = np.nanmin(thetabot[Y<yice])
                cavmax = np.nanmax(thetabot[Y<yice])
                frame = ax1.pcolormesh(X,Y,thetabot,cmap=cmocean.cm.haline,vmin=cavmin,vmax=cavmax)
                cb = plt.colorbar(frame,ax=ax1,pad=0.05)
                cavmin = np.nanmin(thetatop[Y<yice])
                cavmax = np.nanmax(thetatop[Y<yice])
                frame = ax2.pcolormesh(X,Y,thetatop,cmap=cmocean.cm.haline,vmin=cavmin,vmax=cavmax)
                cb = plt.colorbar(frame,ax=ax2,pad=0.05)

    
    ax1.quiver(X[::2,::2],Y[::2,::2],ubot[::2,::2],vbot[::2,::2],scale=0.75)
    ax2.quiver(X[::2,::2],Y[::2,::2],utop[::2,::2],vtop[::2,::2],scale=0.75)

    ax1.set_title("Bottom")
    ax2.set_title("Top")

    if k==1:
        cb = plt.colorbar(frame,ax=ax2,pad=0.05)
    plt.suptitle(description)
    if show:
        plt.show()
    if savepath:
        plt.savefig(savepath)
    plt.close()

def iceFaceMelt(fname,description,times=np.array([])):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    bmask = bottomMask(fname,ds)
    icem = icemask(fname,ds)
    shortname, fpath = outPath(fname) 
    fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)
    theta = ds.THETA.values
    uvel = ds.UVEL.values
    vvel = ds.VVEL.values
    wvel = ds.WVEL.values
    mask = np.full_like(uvel[0,0],1,dtype=bool)
    mask[::4,::4] = 0
    shifwflx = ds.SHIfwFlx.values

    X,Y = np.meshgrid(range(uvel.shape[2]),range(uvel.shape[3]))

    for k in tqdm(range(int(ds.UVEL.values.shape[0]))):
        thetak = theta[k]
        thetabot = thetak *bmask
        thetabot[~bmask] = np.nan
        thetabot = np.nanmean(thetabot,axis=0)
        utop =  np.nanmean(uvel[k] * icem,axis=0)
        vtop =  np.nanmean(vvel[k] * icem,axis=0)
        wtop =  np.nanmean(wvel[k] * icem,axis=0)
        frame = ax1.pcolormesh(utop)
        if k==1:
            cb = plt.colorbar(frame,ax=ax1,pad=0.1)
        #ax1.quiver(ubot,vbot,scale=1.5)

        thetatop = thetak *icem
        thetatop[~icem] = np.nan
        thetatop = np.nanmean(thetatop,axis=0)+1.8
        frame = ax2.pcolormesh(thetatop)
        if k==1:
            cb = plt.colorbar(frame,ax=ax2,pad=0.1)

        ax3.pcolormesh((thetatop**1)*np.sqrt(utop**2+vtop**2+wtop**2))
        ax4.pcolormesh(-shifwflx[k])
        plt.show()
 

def crossSectionAnim(fname,description,quant="THETA",res=1,dim="zonal"):
    fig,ax1 = plt.subplots()
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    extra_variables = dict( KPPdiffS = dict(dims=["k","j","i"], attrs=dict(standard_name="KPPDIFFS", units="kg/m^3")))
    times=getIterNums(fname)
    #ds[quant].values=ds[quant].values*ds.hFacC.values

    #print(ds.hFacC)
    #zonal_average = ds.where(ds.hFacC==1).mean(dim="XC",skipna=True)
    #zonal_average = ds.mean(dim="XC",skipna=True)
    moviewriter = FFMpegFileWriter(fps=1)
    if quant!="DENS":
        ds = open_mdsdataset(fname,prefix=quant,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
        if dim == "zonal":
            #zonal_average = ds.where(ds.hFacC==1).mean(dim="XC",skipna=True)
            zonal_average = ds.isel(XC=95)
            ys = zonal_average.YC.values
        if dim == "meridional":
            zonal_average = ds.isel(YC=76)
            ys = zonal_average.XC.values

        #zonal_average = ds.isel(XC=90)
        zvals = zonal_average[quant].values
        #zvals[zvals==0]=np.nan
        m = np.nanmedian(zvals)
        print(m)
        s = np.nanstd(zvals)
        tmin, tmax = m-5*s,m+s*2
        shortname, fpath = outPath(fname) 
        #plt.hist(zvals[:,:,:].flatten())
        #plt.show()
        fig.suptitle(shortname)
        zs = zonal_average.Z.values
        length = zvals.shape[0]
    if quant=="DENS":
        ds = open_mdsdataset(fname,prefix=["SALT","THETA"],ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
        #zonal_average = ds.where(ds.hFacC==1).isel(XC=100)
        if dim == "zonal":
            zonal_average = ds.where(ds.hFacC==1).mean(dim="XC",skipna=True)
            ys = zonal_average.YC.values
        if dim == "meridional":
            zonal_average = ds.isel(YC=76)
            ys = zonal_average.XC.values
        shortname, fpath = outPath(fname) 
        fig.suptitle(shortname)
        zs = zonal_average.Z.values
        #tmin, tmax = np.nanmin(zonal_average[quant]), np.nanmax(zonal_average[quant])
        zvals = (zonal_average["SALT"].values,zonal_average["THETA"].values)
        length = zvals[0].shape[0]
    with moviewriter.saving(fig, fpath+quant+"|76"+dim+".mp4" , dpi=250):
        print(fpath+quant+"|76"+dim+".mp4")
        for k in tqdm(range(0,length,res)):
            if quant == "DENS":
                frame = ax1.pcolormesh(ys,zs,gsw.sigma0(zvals[0][k,:,:],zvals[1][k,:,:]),cmap="jet",vmin=27.4,vmax=27.7)
            elif quant == "THETA":
                frame = ax1.pcolormesh(ys,zs,zvals[k,:,:],cmap="jet",vmin=-2,vmax=-1.5)
                #frame = ax1.pcolormesh(ys,zs,zvals[k,:,:],cmap="jet",vmin=-2,vmax=1.5)
            elif quant == "SALT":
                frame = ax1.pcolormesh(ys,zs,zvals[k,:,:],cmap="jet",vmin=34,vmax=35.05)
            else:
                frame = ax1.imshow(zvals[k,:,:],cmap="jet",vmin=-0.005,vmax=0.005)
            cb = plt.colorbar(frame)
            moviewriter.grab_frame()
            cb.remove()
            frame.remove()
    plt.close()

def bottomAnim(fname,description,times=np.array([]),quant="THETA",res=5):
    fig,ax1 = plt.subplots()
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    z = np.concatenate((ds.hFacC.values[:-1,:,:]-ds.hFacC.values[1:,:,:],ds.hFacC[-1:,:,:]),axis=0)
    z[z<0]=0
    z[-1,:,:]=0
    zmask = z
    #zonal_average = ds.isel(XC=32)
    moviewriter = FFMpegFileWriter(fps=1)
    tmin, tmax = np.nanmin(ds[quant]), np.nanmax(ds[quant])
    shortname, fpath = outPath(fname) 
    depth = depthFromdZ(ds)
    quantvals = ds[quant].values
    with moviewriter.saving(fig, fpath+"-bot.mp4", dpi=250):
        for k in tqdm([0]+list(range(quantvals.shape[0]))+[-1]):
            d = quantvals[k]
            X = np.full_like(d,np.nan,dtype=float)
            X[ds.hFacC.values != 0]= d[ds.hFacC.values != 0]
            znew = np.multiply(zmask,X)
            nancount = np.nansum(np.isnan(znew),axis=0)
            znew = np.nansum(znew,axis=0)
            znew[nancount==X.shape[0]] = np.nan
            #frame = ax1.pcolormesh(ds.XC.values,ds.YC.values,znew,cmap="jet",vmin=-2,vmax=1)
            frame = ax1.pcolormesh(ds.XC.values,ds.YC.values,znew,cmap="jet",vmin=-2.1,vmax=-1.8)
            CS = ax1.contour(ds.XC.values,ds.YC.values,depth,colors="black",levels=50)
            ax1.clabel(CS, CS.levels, inline=True, fontsize=10)
            cb = plt.colorbar(frame)
            moviewriter.grab_frame()
            cb.remove()
            frame.remove()


def surfaceAnim(fname,description,times=np.array([]),quant="SALT"):
    fig,ax1 = plt.subplots()
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    z = np.concatenate((ds.hFacC.values[:-1,:,:]-ds.hFacC.values[1:,:,:],ds.hFacC[-1:,:,:]),axis=0)
    z[z>0]=0
    z[-1,:,:]=0
    zmask = z
    #zonal_average = ds.isel(XC=32)
    moviewriter = FFMpegFileWriter(fps=1)
    tmin, tmax = np.nanmin(ds[quant]), np.nanmax(ds[quant])
    shortname, fpath = outPath(fname) 
    depth = depthFromdZ(ds)
    quantvals = ds[quant].values#
    with moviewriter.saving(fig, fpath+"-surf.mp4", dpi=250):
        for k in tqdm([0]+list(range(0,quantvals.shape[0]))+[-1]):
            #X = np.full_like(quantvals[k],np.nan,dtype=float)
            #X[ds.hFacC.values != 0]= quantvals[k][ds.hFacC.values != 0]
            #znew = np.multiply(zmask,X)gy
            #nancount = np.nansum(np.isnan(znew),axis=0)
            #znew = np.nansum(znew,axis=0)
            #znew[nancount==X.shape[0]] = np.nan
            znew = quantvals[k][0,:,:]
            frame = ax1.pcolormesh(ds.XC.values,ds.YC.values,znew,cmap="jet",vmin=30)
            ax1.contour(ds.XC.values,ds.YC.values,depth,colors="black",levels=20)
            cb = plt.colorbar(frame)
            moviewriter.grab_frame()
            cb.remove()
            frame.remove()



def folderMap(runsdict,savepath=False):
    fig,axises = plt.subplots(1,1,figsize=(8,7))
    xs,ys,eyeds = [],[],{}
    stats = {"deltaH":[],"Tcdw":[],"gprime":[],"ices":[]}
    statscounter = 0
    for k in runsdict.keys():
        for f in glob.glob(str("../experiments/"+k+"/*"), recursive = True):
            for l in range(len(runsdict[k]["specialstring"])):
                key=runsdict[k]["specialstring"][l]
                if key and "/"+key in f and key == f.rsplit('/', 1)[-1]:
                    try:
                        x,y,newstats=FStheory(f+"/results",None,True)
                        for j in newstats.keys():
                            stats[j].append(newstats[j])
                        if ~np.isnan(y) and ~np.isnan(x):
                            xs.append(x)
                            ys.append(y)
                            eyeds[f+str(l)]=len(xs)-1
                    except:
                        print("yeesh")
                elif not key:
                    try:
                        x,y,newstats=FStheory(f+"/results",None,true)
                        for j in newstats.keys():
                            stats[j].append(newstats[j])
                        if ~np.isnan(y) and ~np.isnan(x):
                            xs.append(x)
                            ys.append(y)
                            eyeds[f+str(l)]=len(xs)-1
                    except:
                        print("yeesh")
    for k in stats.keys():
        print(k,np.nanmean(stats[k]),np.nanstd(stats[k]))
    print(xs,ys)
    xs = np.asarray(([xs])).reshape((-1, 1))
    model = LinearRegression(fit_intercept=False).fit(xs, ys)
    rho0 = 1025
    rhoi = 910
    Cp = 4186
    If = 334000
    C = model.coef_
    print("model_coef: ",model.coef_)
    #W0 = (rho0*Cp)/(rhoi*If*C)
    W0 =  100000#(rho0*Cp)/(rhoi*If*C)
    alpha =  C/((rho0*Cp)/(rhoi*If*W0))
    print("alpha:",alpha)
    oldxs=xs
    #xs= xs*(rho0*Cp)/(rhoi*If*325000)
    xs=model.predict(xs)
    plt.text(.05, .95, '$r^2=$'+str(round(float(pearsonr(oldxs.flatten(),ys)[0])**2,2)), ha='left', va='top', transform=plt.gca().transAxes,fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    plt.gca().set_xlabel(r"$\dot{m}_{\mathrm{pred}} (m/yr)$",fontsize=24)
    plt.gca().set_ylabel(r'$\dot{m}_{\mathrm{obs}} (m/yr)$',fontsize=24)


    plt.plot([0.2,1.4],[0.2,1.4],linestyle="dashed")
    plt.xlim(0,1.4)
    plt.ylim(0.2,1.4)
    for k in runsdict.keys():
        for f in glob.glob(str("../experiments/"+k+"/*"), recursive = True):
            for l in range(len(runsdict[k]["specialstring"])):
                key=runsdict[k]["specialstring"][l]
                if key and key in f:
                    if f+str(l) in eyeds.keys():
                        steadyStateAverageSimple(f+"/results",xs[eyeds[f+str(l)]],fig,axises,color=runsdict[k]["color"][l],marker=runsdict[k]["marker"][l],title=runsdict[k]["description"][0])
                elif not key:
                    try:
                        if f+str(l) in eyeds.keys():
                            steadyStateAverageSimple(f+"/results",xs[eyeds[f+str(l)]],fig,axises,color=runsdict[k]["color"][l],marker=runsdict[k]["marker"][l],title=runsdict[k]["description"][0])
                    except:
                        print("yeesh")
            if savepath:
                plt.savefig(savepath)
    plt.legend()

def folderMapGeneric(func,runsdict,savepath=False,xlabel="",ylabel="",zlabel = "",threed=False):
    if not threed:
        fig,axises = plt.subplots(1,1,figsize=(8,7))
    else:
        fig = plt.gcf()
        axises = plt.figure().add_subplot(projection='3d')

    count = 0
    for k in runsdict.keys():
        for f in glob.glob(str("/jbod/gdf/MITgcm_CS/experiments/"+k+"/*"), recursive = True):
            for l in range(len(runsdict[k]["specialstring"])):
                key=runsdict[k]["specialstring"][l]
                if key and "/"+key in f and key == f.rsplit('/', 1)[-1]:
                    #try:
                    if threed:
                        func(f+"/results",None,fig,axises,color=runsdict[k]["color"][l],marker=runsdict[k]["marker"][l],title=key,znum=count)
                    else:
                        func(f+"/results",None,fig,axises,color=runsdict[k]["color"][l],marker=runsdict[k]["marker"][l],title=key)
                    count+=1
                    #except:
                        #print("yeesh")
                elif not key and False:
                    #try:
                    func(f+"/results",None,fig,axises,color=runsdict[k]["color"][l],marker=runsdict[k]["marker"][l],title=key)
                    #except:
                        #print("yeesh")

            if savepath:
                plt.savefig(savepath)        

    #plt.xlabel(r'Chapman (1998) deep convection depth (m)',fontsize=16)
    #plt.ylabel(r'Average vertical velocity within polyna region across z=250m (m/s)',fontsize=16)

    #plt.xlabel(r'$g^\prime_{\text{disconnected}}s_{\text{ice}} H_{\text{ent}} (\text{Tf}_{\text{z=0}} - \text{Tf}_{\text{z=gl}})$',fontsize=16)
    plt.xlabel(xlabel,fontsize=24)
    plt.ylabel(ylabel,fontsize=24)
    if threed:
        plt.gca().set_zlabel(zlabel,fontsize=24)
    plt.gca().tick_params(axis='both', which='major', labelsize=12)
    #plt.axhline(y=0)
    #plt.axvline(x=-250)
    plt.legend()
    plt.savefig('/jbod/gdf/MITgcm_CS/pics/{}.svg'.format(ylabel))

def folderMapMoreGeneric(func,runsdict):
    for k in runsdict.keys():
        for f in glob.glob(str("/home/garrett/Projects/MITgcm_CS/experiments/"+k+"/*"), recursive = True):
            for l in range(len(runsdict[k]["specialstring"])):
                print(f)
                key=runsdict[k]["specialstring"][l]
                if "at125" in f:
                    if key and key in f:
                        #try:
                        func(f+"/results",key)
                        #except:
                        #print("yeesh")
                    elif not key:
                        #try:
                        func(f+"/results",key)
                        #except:
                            #print("yeesh")

def folderMapRefresh(runsdict,save=False):
    prepath = os.path.abspath(os.getcwd()).replace("analysis","experiments")
    for k in runsdict.keys():
        for f in tqdm(glob.glob(str(prepath+"/"+k+"/*"), recursive = True)):
            for l in range(len(runsdict[k]["specialstring"])):
                key=runsdict[k]["specialstring"][l]
                if key and key == f.rsplit('/', 1)[-1] :
                    #try:
                    timeSeries(f+"/results",True)
                    print(f)
                    #except:
                        #print("yeesh")
                #elif not key:
                    #try:
                        #if f+str(l) in eyeds.keys():
                            #timeSeries(f+"/results",True)
                    #except:
                        #print("yeesh")
#
            #if save:
                #plt.savefig("/home/garrett/Projects/HUB/paperfigures/"+k+".png")
def meltMapAverage(fname,description,res=1,ax1=None,show=False,savepath=False):
    if not ax1:
        fig,ax1 = plt.subplots(figsize=(10,8))
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,prefix=["SHIfwFlx"],ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    mask = np.logical_and(ds.hFacC.values[0]==0,np.sum(ds.hFacC.values,axis=0)!=0)
    times=np.asarray(times)
    times = times*grabDeltaT(fname)/60.0/60.0/24.0/365.0
    melt= -np.nanmean(ds.SHIfwFlx.values[times>2.5],axis=0)*(60*60*24*365)*(1/920.0)
    xs = ds.XG.values/1000
    ys = ds.YG.values/1000
    newcmap = cmocean.tools.crop(cmocean.cm.balance, 0, 3, 0)

    vals = grabMatVars(fname,("h","icedraft"))
    icedraft = np.abs(np.asarray(vals["icedraft"]))
    h = np.abs(np.asarray(vals["h"]))
    icedraft = np.logical_and(icedraft!=0,icedraft!=h)
    melt[~icedraft.T]=np.nan
    
    im = ax1.pcolormesh(xs,ys,melt,cmap=cmocean.cm.balance,vmin=-5,vmax=5)

    bound = np.argwhere(~np.isnan(melt))

    #ax2.set_xlim((140,260))
    #ax2.set_ylim((0,160))
    pad = 10
    ax1.set_xlim(xs[min(bound[:, 1])]-pad, xs[max(bound[:, 1])]+pad)
    ax1.set_ylim(ys[min(bound[:, 0])]-pad, ys[max(bound[:, 0])]+pad)

    ax1.set_xlabel("x (km)",fontsize=18)
    ax1.set_ylabel("y (km)",fontsize=18)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    caxout = inset_axes(
        ax1,
        width="2%",  # width: 5% of parent_bbox width
        height="40%",  # height: 50%
        loc="lower left",
        bbox_to_anchor=(1.0, 0.30, 1, 1),
        bbox_transform=ax1.transAxes,
        borderpad=1,
    )
    caxout.tick_params(labelsize=18)
    plt.colorbar(im,cax=caxout,ticks=[0,3])
    caxout.set_ylabel('m/yr', rotation=0,fontsize=18)

    ax1.set_ylabel('Y (km)',fontsize=18)
    ax1.set_xlabel('X (km)',fontsize=18)
    plt.title(description)
    if show:
        plt.show() 
    if savepath:
        plt.savefig(savepath)
    #caxout = plt.colorbar(im,ax=ax1, aspect=40,shrink=0.4,ticks=range(0,41,10))
    #caxout.ax.tick_params(labelsize=18)
def crossAndMelt(fname,name=""):
    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(15,6))
    plt.subplots_adjust(wspace=0.45)
    crossSectionAverage(fname,"Reference",quant="THETA",dim="zonal",ax1=ax1,show=False)
    meltMapAverage(fname,"",ax1=ax2)
    plt.show()
    plt.savefig("/home/garrett/Projects/HUB/paperfigures/crossAndMelts/"+name+".png")
    plt.close()

def folderMapTimeSeries(runsdict,save=True):
    fig,axises = plt.subplots(2,3,figsize=(8,7))
    for k in runsdict.keys():
        for f in glob.glob(str("../experiments/"+k+"/*"), recursive = True):
            for l in range(len(runsdict[k]["specialstring"])):
                key=runsdict[k]["specialstring"][l]
                if key and key == f.rsplit('/', 1)[-1]:
                    #try:
                        #print(key)
                    print(f)
                    timeSeriesDashboard(f+"/results",key+f[-6:-10],fig,axises,color=runsdict[k]["color"][l])
                    #except:
                        #print("yeesh")
    axises[0][0].legend()

    plt.show()


def TSAnim(fname,description,res=1):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    extra_variables = dict( KPPdiffS = dict(dims=["k","j","i"], attrs=dict(standard_name="KPPDIFFS", units="kg/m^3")))
    times=getIterNums(fname)
    #ds[quant].values=ds[quant].values*ds.hFacC.values

    #print(ds.hFacC)
    #zonal_average = ds.where(ds.hFacC==1).mean(dim="XC",skipna=True)
    #zonal_average = ds.mean(dim="XC",skipna=True)
    moviewriter = FFMpegFileWriter(fps=1)
    ds = open_mdsdataset(fname,prefix=["SALT","THETA"],ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    salt = ds.SALT.values
    theta = ds.THETA.values
    length = salt.shape[0]
    shortname, fpath = outPath(fname) 
    fig,ax1 = plt.subplots()
    with moviewriter.saving(fig, fpath+"TS"+".mp4" , dpi=250):
        print(fpath+"TS"+".mp4")
        for k in tqdm(range(0,length,res)):
            frame = ax1.scatter(salt[k].flatten()[mask],theta[k].flatten()[mask],s=0.5,color="black")
            ax1.set_xlim([34.1,34.75])
            ax1.set_ylim([-2.25,-1.790])
            cb = plt.colorbar(frame)
            moviewriter.grab_frame()
            cb.remove()
            frame.remove()
    plt.close()


def volumetricTS(fname,description,res=1,show=True,savepath=False):
    data = timeSeries(fname)
    for k in data.keys():
        if k != "ts":
            try:
                data[k] = np.nanmean(data[k][data["ts"]>2])
            except:
                data[k]=np.nan
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    extra_variables = dict( KPPdiffS = dict(dims=["k","j","i"], attrs=dict(standard_name="KPPDIFFS", units="kg/m^3")))
    times=getIterNums(fname)
    #ds[quant].values=ds[quant].values*ds.hFacC.values

    #print(ds.hFacC)
    #zonal_average = ds.where(ds.hFacC==1).mean(dim="XC",skipna=True)
    #zonal_average = ds.mean(dim="XC",skipna=True)
    ds = open_mdsdataset(fname,prefix=["SALT","THETA"],ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    yice = grabMatVars(fname[:-1],("Yicefront"))["Yicefront"][0][0]
    variables = grabMatVars(fname,('tNorth','tEast','sEast','zz','pp'))

    tEast = np.asarray(variables["tEast"])
    sEast = np.asarray(variables["sEast"])
    tEast = tEast[int(tEast.shape[0]-1)]
    sEast = sEast[int(sEast.shape[0]-1)]

    # ds = ds.where(ds.YC<yice-10000,drop=True)
    ds = ds.where(ds.YC<150000,drop=True)
    # ds = ds.where(ds.YC>yice-5000,drop=True)
    # ds = ds.where(ds.XC>150000,drop=True)
    print((yice-40000)/1000)
    salt = ds.SALT.values
    theta = ds.THETA.values
    for i in range(salt.shape[0]):
        salt[i,ds.hFacC.values<0.1] = np.nan
        theta[i,ds.hFacC.values<0.1] = np.nan
    Z,Y,X = np.meshgrid(ds.Z.values,ds.YC.values,ds.XC.values,indexing='ij')
    volume = (ds.hFacC * ds.drF * ds.rA).values
    length = salt.shape[0]
    shortname, fpath = outPath(fname) 
    fig,ax1 = plt.subplots(figsize=(10,8))

    times=np.asarray(times)
    times = times*grabDeltaT(fname)/60.0/60.0/24.0/365.0
    print(times)
    timemask = times>8


    Spolyna = data["scdw"]#np.max(np.mean(salt[:10],axis=0)[volume>0])
    Tpolyna = data['tcdw']#np.max(np.mean(theta[:10],axis=0)[volume>0])
    from analysis import fpAtGl
    meltdepth = 0
    Tf = fpAtGl(meltdepth,Spolyna)
    import gsw
    print(gsw.pt_from_t(Spolyna,Tf,meltdepth,0),Tf)
    ax1.axhline(y=gsw.pt_from_t(Spolyna,Tf,meltdepth,0))
    Cp = 3390
    If = 334000

    Tf = gsw.pt_from_t(Spolyna,Tf,meltdepth,0)
    Sm = Spolyna/(1-(Cp/If)*(Tf-(Tpolyna)))

    # im=ax1.hist2d(np.mean(salt[timemask],axis=0).flatten(),np.mean(theta[timemask],axis=0).flatten(),weights=volume.flatten(),density=True,bins=100,cmin=0.001,range=[[34,35],[-2.4,-1.6]])#,norm="log")
    meansalt = np.mean(salt[timemask],axis=0)
    meantemp = np.mean(theta[timemask],axis=0)
    tfs = fpAtGl(Z,meansalt)
    tfs = gsw.pt_from_t(meansalt,tfs,Z,0)
    #im=ax1.scatter(meansalt,np.mean(theta[timemask],axis=0).flatten(),c=Z.flatten(),vmin=-900,vmax=-200)#,norm="log")
    im=ax1.scatter(meansalt,np.mean(theta[timemask],axis=0).flatten(),c=Z.flatten(),vmin=-900,vmax=0)#,norm="log")
    # im=ax1.scatter(meansalt.flatten(),meantemp.flatten(),c=tfs.flatten()-meantemp.flatten()) #,norm="log")
    plt.colorbar(im,ax=ax1)

    frac = (data['ssurf']-Sm)/(Spolyna-Sm)
    print(frac)

    Cp = 3390
    If = 334000

    for meltdepth in range(100,1000,100):
        Tf = fpAtGl(meltdepth,Spolyna)
        Tf = gsw.pt_from_t(Spolyna,Tf,meltdepth,0)
        print(Tf)
        print((1-(Cp/If)*(Tf-(Tpolyna))))
        Sm = Spolyna/(1-(Cp/If)*(Tf-(Tpolyna)))
        Smid = Sm
        Tmid = Tf
        ax1.plot((Smid,Spolyna),(Tmid,Tpolyna),c="red")


    # plt.scatter(((Smid+Spolyna)/2,(Sm+Sm)/2),((Tmid+Tpolyna)/2,(Tf+Tf)/2),c="purple")
    print(data["scdw"]-data["ssurf"])
    plt.scatter((data["ssurf"],data["scdw"]),(data["tsurf"],data['tcdw']),c="red")
    #plt.colorbar(im[3],ax=ax1)
    #ax1.set_xlim([34.1,34.75])
    #ax1.set_ylim([-2.25,-1.790])
    plt.title(description)

    if show:
        plt.show()
    if savepath:
        plt.savefig(savepath)

def buildPortfolio(fname,name):
    shortname, fpath = outPath(fname) 
    foliopath = "../pics/"+name+"-portfolio"
    if not os.path.exists(foliopath):
        os.makedirs(foliopath)
    volumetricTS(fname,shortname,show=False,savepath=foliopath+"/volumetricTS.png")
    meltMapAverage(fname,shortname,show=False,savepath=foliopath+"/meltmap.png")
    crossSectionAverage(fname,shortname,150*10**3,quant="THETA",dim="meridional",show=False,savepath=foliopath+"/150meridionalT.png")
    crossSectionAverage(fname,shortname,150*10**3,quant="THETA",dim="meridional",show=False,savepath=foliopath+"/150meridionalT-nf.png",fixcb=False)
    crossSectionAverage(fname,shortname,200*10**3,quant="THETA",dim="meridional",show=False,savepath=foliopath+"/200meridionalT.png")
    crossSectionAverage(fname,shortname,200*10**3,quant="THETA",dim="meridional",show=False,savepath=foliopath+"/200meridionalT-nf.png",fixcb=False)
    crossSectionAverage(fname,shortname,250*10**3,quant="THETA",dim="meridional",show=False,savepath=foliopath+"/250meridionalT.png")
    crossSectionAverage(fname,shortname,250*10**3,quant="THETA",dim="meridional",show=False,savepath=foliopath+"/250meridionalT-nf.png",fixcb=False)
###
    crossSectionAverage(fname,shortname,150*10**3,quant="SALT",dim="meridional",show=False,savepath=foliopath+"/150meridionalS.png")
    crossSectionAverage(fname,shortname,150*10**3,quant="SALT",dim="meridional",show=False,savepath=foliopath+"/150meridionalS-nf.png",fixcb=False)
    crossSectionAverage(fname,shortname,200*10**3,quant="SALT",dim="meridional",show=False,savepath=foliopath+"/200meridionalS.png")
    crossSectionAverage(fname,shortname,200*10**3,quant="SALT",dim="meridional",show=False,savepath=foliopath+"/200meridionalS-nf.png",fixcb=False)
    crossSectionAverage(fname,shortname,250*10**3,quant="SALT",dim="meridional",show=False,savepath=foliopath+"/250meridionalS.png")
    crossSectionAverage(fname,shortname,250*10**3,quant="SALT",dim="meridional",show=False,savepath=foliopath+"/250meridionalS-nf.png",fixcb=False)
#

    crossSectionAverage(fname,shortname,150*10**3,quant="UVEL",dim="meridional",show=False,savepath=foliopath+"/150meridionalU.png")
    crossSectionAverage(fname,shortname,150*10**3,quant="UVEL",dim="meridional",show=False,savepath=foliopath+"/150meridionalU-nf.png",fixcb=False)
    crossSectionAverage(fname,shortname,200*10**3,quant="UVEL",dim="meridional",show=False,savepath=foliopath+"/200meridionalU.png")
    crossSectionAverage(fname,shortname,200*10**3,quant="UVEL",dim="meridional",show=False,savepath=foliopath+"/200meridionalU-nf.png",fixcb=False)
    crossSectionAverage(fname,shortname,250*10**3,quant="UVEL",dim="meridional",show=False,savepath=foliopath+"/250meridionalU.png")
    crossSectionAverage(fname,shortname,250*10**3,quant="UVEL",dim="meridional",show=False,savepath=foliopath+"/250meridionalU-nf.png",fixcb=False)


    crossSectionAverage(fname,shortname,150*10**3,quant="DENS",dim="meridional",show=False,savepath=foliopath+"/150meridionalD.png")
    crossSectionAverage(fname,shortname,150*10**3,quant="DENS",dim="meridional",show=False,savepath=foliopath+"/150meridionalD-nf.png",fixcb=False)
    crossSectionAverage(fname,shortname,200*10**3,quant="DENS",dim="meridional",show=False,savepath=foliopath+"/200meridionalD.png")
    crossSectionAverage(fname,shortname,200*10**3,quant="DENS",dim="meridional",show=False,savepath=foliopath+"/200meridionalD-nf.png",fixcb=False)
    crossSectionAverage(fname,shortname,250*10**3,quant="DENS",dim="meridional",show=False,savepath=foliopath+"/250meridionalD.png")
    crossSectionAverage(fname,shortname,250*10**3,quant="DENS",dim="meridional",show=False,savepath=foliopath+"/250meridionalD-nf.png",fixcb=False)
##
    crossSectionAverage(fname,shortname,80*10**3,quant="THETA",dim="zonal",show=False,savepath=foliopath+"/80zonalT.png")
    crossSectionAverage(fname,shortname,80*10**3,quant="THETA",dim="zonal",show=False,savepath=foliopath+"/80zonalT-nf.png",fixcb=False)
    crossSectionAverage(fname,shortname,140*10**3,quant="THETA",dim="zonal",show=False,savepath=foliopath+"/140zonalT.png")
    crossSectionAverage(fname,shortname,140*10**3,quant="THETA",dim="zonal",show=False,savepath=foliopath+"/140zonalT-nf.png",fixcb=False)
    crossSectionAverage(fname,shortname,175*10**3,quant="THETA",dim="zonal",show=False,savepath=foliopath+"/175zonalT.png")
    crossSectionAverage(fname,shortname,175*10**3,quant="THETA",dim="zonal",show=False,savepath=foliopath+"/175zonalT-nf.png",fixcb=False)
##
    crossSectionAverage(fname,shortname,80*10**3,quant="SALT",dim="zonal",show=False,savepath=foliopath+"/80zonalS.png")
    crossSectionAverage(fname,shortname,80*10**3,quant="SALT",dim="zonal",show=False,savepath=foliopath+"/80zonalS-nf.png",fixcb=False)
    crossSectionAverage(fname,shortname,140*10**3,quant="SALT",dim="zonal",show=False,savepath=foliopath+"/140zonalS.png")
    crossSectionAverage(fname,shortname,140*10**3,quant="SALT",dim="zonal",show=False,savepath=foliopath+"/140zonalS-nf.png",fixcb=False)
    crossSectionAverage(fname,shortname,175*10**3,quant="SALT",dim="zonal",show=False,savepath=foliopath+"/175zonalS.png")
    crossSectionAverage(fname,shortname,175*10**3,quant="SALT",dim="zonal",show=False,savepath=foliopath+"/175zonalS-nf.png",fixcb=False)
##
    crossSectionAverage(fname,shortname,80*10**3,quant="DENS",dim="zonal",show=False,savepath=foliopath+"/80zonalD.png")
    crossSectionAverage(fname,shortname,80*10**3,quant="DENS",dim="zonal",show=False,savepath=foliopath+"/80zonalD-nf.png",fixcb=False)
    crossSectionAverage(fname,shortname,140*10**3,quant="DENS",dim="zonal",show=False,savepath=foliopath+"/140zonalD.png")
    crossSectionAverage(fname,shortname,140*10**3,quant="DENS",dim="zonal",show=False,savepath=foliopath+"/140zonalD-nf.png",fixcb=False)
    crossSectionAverage(fname,shortname,175*10**3,quant="DENS",dim="zonal",show=False,savepath=foliopath+"/175zonalD.png")
    crossSectionAverage(fname,shortname,175*10**3,quant="DENS",dim="zonal",show=False,savepath=foliopath+"/175zonalD-nf.png",fixcb=False)

    topMap(fname,shortname,quant="THETA",show=False,savepath=foliopath+"/topbotT.png")
    topMap(fname,shortname,quant="THETA",show=False,savepath=foliopath+"/topbotT-nf.png",fixcb=False)
    topMap(fname,shortname,quant="SALT",show=False,savepath=foliopath+"/topbotS.png")
    topMap(fname,shortname,quant="SALT",show=False,savepath=foliopath+"/topbotS-nf.png",fixcb=False)

def overturning_plot(fname,name):

    times=getIterNums(fname)

    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    ds = open_mdsdataset(fname,prefix=["VVEL"],ignore_unknown_vars=True,iters=times,extra_variables = extra_variables)
    times = np.asarray( [0]+list(ds.time.values))*(10**-9)

    yice = grabMatVars(fname,("Yicefront"))["Yicefront"][0][0]
    yshelf = grabMatVars(fname,("Yshelfbreak"))["Yshelfbreak"][0][0]
    ycoast = grabMatVars(fname,("Ycoast"))["Ycoast"][0][0]
    xeast = grabMatVars(fname,("Xeast"))["Xeast"][0][0]
    xwest = grabMatVars(fname,("Xwest"))["Xwest"][0][0]
    saltfluxvals = grabMatVars(fname,("saltfluxvals"))["saltfluxvals"].T

    ices = slope(fname)

    vvel = ds.VVEL

    transport = (vvel*ds.hFacS*ds.drF*ds.dxG)[10:].mean(dim="time").sum(dim="XC").where(ds.YG<190000)
    

    transport = transport.cumsum(dim="Z")
    transport = transport.where(ds.hFacS.sum(dim="XC")>0)

    streamlow = transport.where(np.abs(transport.YG-(yice+15000))<6000).min().values
    streamhigh = transport.where(np.abs(transport.YG-(yice+15000))<6000).max().values
    print(streamlow)
    
    transport['YG'] = transport.YG/1000
    contours = transport.plot.contour(levels=np.linspace(streamlow+100,streamhigh-100,num = 10 ),colors="gray",alpha=1,aspect=1.5,size=8)
    im = (transport/(1e6)).plot.pcolormesh(label="$Stream function (Sv)$",add_colorbar=False,cmap="cmo.balance")
    start_xs = range((yice+15000)-6000,(yice+15000)-6000,10)
    start_ys = range(-590,0,10)
    start_X,start_Y = np.meshgrid(start_xs,start_ys)
    start_points = np.vstack([start_X.ravel(), start_Y.ravel()])
    X,Y = np.meshgrid(transport.YG,np.abs(transport.Z))
    hfacc = ds.hFacC.sum(axis=2).values
    mask = ~(hfacc !=0)
    mask = mask.astype(float)
    mask[mask==0]=np.nan
    mask[np.logical_and((-X*ices*1000 + 900)<Y,~np.isnan(mask))] = 2
    mask[np.logical_and((-X*ices*1000 + 900)>Y,~np.isnan(mask))] = -2
    cmap = ListedColormap(["#e4f3f8", "#cab590"])

    plt.pcolormesh(transport.YG,transport.Z,mask,cmap=cmap)
    # pdb.set_trace()
    # plt.gca().streamplot(X,Y,-np.diff(transport.values,axis=1)/np.diff(Y),-np.diff(transport.values,axis=0,prepend=0)/np.diff(X,prepend=0),start_points=start_points)
    plt.xlim(0,190)
    plt.ylim(-1000,0)
    plt.gca().set_xlabel('y (km)', fontsize=24)
    plt.gca().set_ylabel('z (m)', fontsize=24)
    plt.gca().tick_params(axis='both', labelsize=16)
    cbar = plt.colorbar(im,label="$Stream function (Sv)$")
    cbar.ax.yaxis.label.set_size(24)
    cbar.ax.tick_params(axis='both', labelsize=16)

    plt.tight_layout()
    plt.savefig('/jbod/gdf/MITgcm_CS/pics/overturning'+name+".png",dpi=350)

    # bools = []
    # for i in range(len(contours.collections)):
    #     if len(contours.collections[i].get_paths())>0:
    #         xmax = contours.collections[i].get_paths()[0].vertices[:,0].max()
    #         xmin = contours.collections[i].get_paths()[0].vertices[:,0].min()
    #         bools.append((xmax>yice+20000) and (xmin<yice))
    
    # bools = np.asarray(bools)
    # print(bools)
    # print(np.sum(bools)/len(bools))
    plt.show()

