#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 18:44:59 2017

@author: pramos

Usufull functions for UVW velocity-space treatments
"""
import numpy as np


#constants
incl=np.deg2rad(62.87124882)    #inclination of galactic plane
alom=np.deg2rad(282.8594813)    #RA of the equatorial node
lom=np.deg2rad(32.93680516)     #longitude gal of Ascending node of the galactic equator
deo=np.deg2rad(27.12825111)     #dec NGP 
alpm=np.deg2rad(192.8594813)    #RA NGP
theta=np.deg2rad(122.93680516)  #londitude NCP
const=4.7404705                 #conversion factor km·yr/s

def radec2lb(ra,dec):
    """
    Returns galactic longitude and latitude (radians) given
    the right ascention and declination (radians)
    """
    #latitude
    sinb=np.sin(dec)*np.cos(incl)-np.cos(dec)*np.sin(incl)*np.sin(ra-alom)
    b=np.arcsin(sinb)
    #longitude
    sinl=(np.sin(dec)*np.sin(incl)+np.cos(dec)*np.cos(incl)*np.sin(ra-alom))/np.cos(dec)
    cosl=np.cos(dec)*np.cos(ra-alom)/np.cos(dec)
    l=np.arctan2(sinl,cosl)+lom
    if l<0:l=l+2*np.pi
    
    return l,b

def pmradec2lb(ra,dec,l,b,mua,mud):
    """
    Returns proper motion in galactic coordinates given
    the right ascention, declination,Gal.Lat and Long. (radians) and
    equatorial proper motions (pmra,pmdec in mas/yr)
    """
    siphi=np.cos(deo)*np.sin(ra-alpm)/np.cos(b)
    cophi=(np.sin(deo)-np.sin(dec)*np.sin(b))/(np.cos(dec)*np.cos(b))
    
    return mua*cophi+mud*siphi,-mua*siphi+mud*cophi

def Jacob(x):
    """
    Returns the Jacobian of the transformation from ICRS to l,b,plx,U,V,W.
    The input parameters is an iterator of the form:
        x[0]: right ascention (equatorial) -> degrees
        x[1]: declination (equatorial) -> degrees
        x[2]: parallax -> mas
        x[3]: proper motion (mualpha*) -> mas/yr
        x[4]: proper motion (mudelta) -> mas/yr
        x[5]: radial velocity -> km/s
    """
    
    #inicialization
    jac0=np.zeros([6,6])
    jac1=np.zeros([6,6])
    jac2=np.zeros([6,6])
    jac3=np.zeros([6,6])
    #RA & DEC in degrees
    a=np.deg2rad(x[0])
    d=np.deg2rad(x[1])
    #Parallax in mas
    w=x[2]
    #proper motion in mas/yr
    mua=x[3]
    mud=x[4]
    #radial velocity along l.o.s. in km/s
    vrad=x[5]
    #galactic coordinates in degrees
    l,b=radec2lb(a,d)

    cd=np.cos(d);sd=np.sin(d);tand=np.tan(d)
    cb=np.cos(b); sb=np.sin(b);tanb=np.tan(b)
    cl=np.cos(l); sl=np.sin(l)
    cdeo=np.cos(deo); sdeo=np.sin(deo)
    cincl=np.cos(incl); sincl=np.sin(incl)
    calom=np.cos(a-alom); salom=np.sin(a-alom); talom=np.tan(a-alom)
    clt=np.cos(l-theta); slt=np.sin(l-theta)
    
    #units
    jac0[0,0]=np.pi/180/3600
    jac0[1,1]=np.pi/180/3600
    jac0[2,2]=1
    jac0[3,3]=1
    jac0[4,4]=1
    jac0[5,5]=1
    
    #from equatorial to galactic coordinates
    
    #l-alpha
    jac1[0,0]=((cincl/(calom)**2+sincl/calom*talom*tand)/(1+(cincl*talom+sincl/calom*tand)**2))
    #l-delta
    jac1[0,1]=(sincl/(calom*(cd)**2)/(1+(cincl*talom+sincl/calom*tand)**2))
    #b-alpha
    jac1[1,0]=(-((calom*cd*sincl)/np.sqrt(1-(cincl*sd-cd*salom*sincl)**2)))
    
    #b-delta
    jac1[1,1]=((cd*cincl+salom*sd*sincl)/np.sqrt(1-(cincl*sd-cd*salom*sincl)**2))
        
    #parallax-parallax
    jac1[2,2]=1
    
    #mua - mua
    jac1[3,3]=1
    
    #mud - mud
    jac1[4,4]=1
    
    #vrad - vrad
    jac1[5,5]=1
    
    """"""""""""""""""
    #from mua/mub to mul/mub
    
    #l-l
    jac2[0,0]=1
    #b-b
    jac2[1,1]=1
    #par-par
    jac2[2,2]=1
    
    #mul-l
    jac2[3,0]=-((mud*cdeo*clt)/(1-(cb*cdeo*clt+sb*sdeo)**2)**(0.5))-(mua*cdeo*( \
        cb*cdeo*clt+sb*sdeo)*(sdeo-sb*(cb*cdeo*clt+sb*sdeo))*slt)/(1-(cb*cdeo*clt+sb*sdeo)**2)**( \
        1.5)+(mua*cdeo*sb*slt)/(1-(cb*cdeo*clt+sb*sdeo)**2)**(0.5)+(mud*cb*cdeo**2*(cb*cdeo*clt+sb*sdeo)*slt*slt)/( \
    1-(cb*cdeo*clt+sb*sdeo)**2)**(1.5)
    
    #ml-b
    jac2[3,1]=(mua/(cb)*(-(cdeo*clt*sb)+cb*sdeo)*(cb*cdeo*clt+sb*sdeo)*(sdeo-sb*(cb*cdeo*clt+sb*sdeo)))/( \
        1-(cb*cdeo*clt+sb*sdeo)**2)**(1.5)+(mua/(cb)*(-(sb*(-(cdeo*clt*sb)+cb*sdeo))-cb*(cb*cdeo*clt+sb*sdeo)))/ \
    (1-(cb*cdeo*clt+sb*sdeo)**2)**(0.5)-(mud*cdeo*(-(cdeo*clt*sb)+cb*sdeo)*(cb*cdeo*clt+sb*sdeo)*slt)/( \
        1-(cb*cdeo*clt+sb*sdeo)**2)**(1.5)+(mua/(cb)*(sdeo-sb*(cb*cdeo*clt+sb*sdeo))*tanb)/(1-(cb*cdeo*clt+sb*sdeo)**2)**(0.5)
    
    #ml-mua
    jac2[3,3]=(1/(cb)*(sdeo-sb*(cb*cdeo*clt+sb*sdeo)))/(1-(cb*cdeo*clt+sb*sdeo)**2)**(0.5)
    
    #ml-mud
    jac2[3,4]=-((cdeo*slt)/(1-(cb*cdeo*clt+sb*sdeo)**2)**(0.5))
    
    
    #mb-l
    jac2[4,0]=(mua*cdeo*clt)/(1-(cb*cdeo*clt+sb*sdeo)**2)**(0.5)-(mud*cdeo*(cb*cdeo*clt+sb*sdeo)*(sdeo-sb*( \
        cb*cdeo*clt+sb*sdeo))*slt)/(1-(cb*cdeo*clt+sb*sdeo)**2)**(1.5)+(mud*cdeo*sb*slt)/(1-(cb*cdeo*clt+sb*sdeo)**2)**(0.5)-( \
        mua*cb*cdeo**2*(cb*cdeo*clt+sb*sdeo)*slt*slt)/(1-(cb*cdeo*clt+sb*sdeo)**2)**(1.5)
    
    #mb-b
    jac2[4,1]=(mud/(cb)*(-(cdeo*clt*sb)+cb*sdeo)*(cb*cdeo*clt+sb*sdeo)*(sdeo-sb*(cb*cdeo*clt+sb*sdeo)))/ \
    (1-(cb*cdeo*clt+sb*sdeo)**2)**(1.5)+(mud/(cb)*(-(sb*(-(cdeo*clt*sb)+cb*sdeo))-cb*(cb*cdeo*clt+sb*sdeo)))/(1-( \
        cb*cdeo*clt+sb*sdeo)**2)**(0.5)+(mua*cdeo*(-(cdeo*clt*sb)+cb*sdeo)*(cb*cdeo*clt+sb*sdeo)*slt)/( \
        1-(cb*cdeo*clt+sb*sdeo)**2)**(1.5)+(mud/(cb)*(sdeo-sb*(cb*cdeo*clt+sb*sdeo))*tanb)/(1-(cb*cdeo*clt+sb*sdeo)**2)**(0.5)
    
    #mb-mua
    jac2[4,3]=(cdeo*slt)/(1-(cb*cdeo*clt+sb*sdeo)**2)**(0.5)
         
    #mb-mud
    jac2[4,4]=(1/(cb)*(sdeo-sb*(cb*cdeo*clt+sb*sdeo)))/(1-(cb*cdeo*clt+sb*sdeo)**2)**(0.5)
    
    #vrad-vrad
    jac2[5,5]=1
    
    
    "calculating ml and mb"
    
    ml,mb=pmradec2lb(a,d,l,b,mua,mud)
    
    #ml/mb -> uvw
    
    #l-l
    jac3[0,0]=1
    #b-b
    jac3[1,1]=1
    #par-par
    jac3[2,2]=1
    #U-l
    jac3[3,0]=-((const*ml*cl)/w)-vrad*cb*sl+(const*mb*sb*sl)/w
    #U-b
    jac3[3,1]=-((const*mb*cb*cl)/w)-vrad*cl*sb
    
    #U-par
    jac3[3,2]=(const*mb*cl*sb)/w**2+(const*ml*sl)/w**2
    
    #U-ml
    jac3[3,3]=-((const*sl)/w)
    
    #U-mb
    jac3[3,4]=-((const*cl*sb)/w)
    
    #U-vrad
    jac3[3,5]=cb*cl
    
    
    #V-l
    jac3[4,0]=vrad*cb*cl-(const*mb*cl*sb)/w-(const*ml*sl)/w
         
    #V-b
    jac3[4,1]=-((const*mb*cb*sl)/w)-vrad*sb*sl
    
    #V-par
    jac3[4,2]=-((const*ml*cl)/w**2)+(const*mb*sb*sl)/w**2
    
    #V-ml
    jac3[4,3]=(const*cl)/w
    
    #V-mb
    jac3[4,4]=-((const*sb*sl)/w)
    
    #V-vrad
    jac3[4,5]=cb*sl
    
    
    #W-l
    jac3[5,0]=0
    #W-b
    jac3[5,1]=vrad*cb-(const*mb*sb)/w
    #W-par
    jac3[5,2]=-((const*mb*cb)/w**2)
    #W-ml
    jac3[5,3]=0
    #W-mb
    jac3[5,4]=(const*cb)/w
    #W-vrad
    jac3[5,5]=sb
    
    
    return np.dot(jac3,np.dot(jac2,np.dot(jac1,jac0)))


def Jacob4(x):
    """
    Returns the Jacobian (4x4) of the transformation from ICRS to l,b,plx,U,V,W, ignoring positional errors
    to increase speed.
    The input parameters is an iterator of the form:
        x[0]: right ascention (equatorial) -> degrees
        x[1]: declination (equatorial) -> degrees
        x[2]: parallax -> mas
        x[3]: proper motion (mualpha*) -> mas/yr
        x[4]: proper motion (mudelta) -> mas/yr
        x[5]: radial velocity -> km/s
    """
    
    #inicialization
    jac0=np.zeros([4,4])
    #jac1=np.zeros([4,4])
    jac2=np.zeros([4,4])
    jac3=np.zeros([4,4])
    #RA & DEC in degrees
    a=np.deg2rad(x[0])
    d=np.deg2rad(x[1])
    #Parallax in mas
    w=x[2]
    #proper motion in mas/yr
    mua=x[3]
    mud=x[4]
    #radial velocity along l.o.s. in km/s
    vrad=x[5]
    #galactic coordinates in degrees
    l,b=radec2lb(a,d)

    cd=np.cos(d);sd=np.sin(d);tand=np.tan(d)
    cb=np.cos(b); sb=np.sin(b);tanb=np.tan(b)
    cl=np.cos(l); sl=np.sin(l)
    cdeo=np.cos(deo); sdeo=np.sin(deo)
    cincl=np.cos(incl); sincl=np.sin(incl)
    calom=np.cos(a-alom); salom=np.sin(a-alom); talom=np.tan(a-alom)
    clt=np.cos(l-theta); slt=np.sin(l-theta)
    
    #mas -> rad
    jac0[0,0]=1 #par
    jac0[1,1]=1 #pmra
    jac0[2,2]=1 #pmdec
    jac0[3,3]=1 #vr
    
    #from mua/mub to mul/mub
    
    #par-par
    jac2[0,0]=1
    
    #ml-mua
    jac2[1,1]=(1/(cb)*(sdeo-sb*(cb*cdeo*clt+sb*sdeo)))/(1-(cb*cdeo*clt+sb*sdeo)**2)**(0.5)
    
    #ml-mud
    jac2[1,2]=-((cdeo*slt)/(1-(cb*cdeo*clt+sb*sdeo)**2)**(0.5))
    
    
    #mb-mua
    jac2[2,1]=(cdeo*slt)/(1-(cb*cdeo*clt+sb*sdeo)**2)**(0.5)
         
    #mb-mud
    jac2[2,2]=(1/(cb)*(sdeo-sb*(cb*cdeo*clt+sb*sdeo)))/(1-(cb*cdeo*clt+sb*sdeo)**2)**(0.5)
    
    #vrad-vrad
    jac2[3,3]=1
    
    
    "calculating ml and mb"
    
    ml,mb=pmradec2lb(a,d,l,b,mua,mud)
    
    #ml/mb -> uvw
    
    #par-par
    jac3[0,0]=1
    
    #U-par
    jac3[1,0]=(const*mb*cl*sb)/w**2+(const*ml*sl)/w**2
    
    #U-ml
    jac3[1,1]=-((const*sl)/w)
    
    #U-mb
    jac3[1,2]=-((const*cl*sb)/w)
    
    #U-vrad
    jac3[1,3]=cb*cl
    
    
    #V-par
    jac3[2,0]=-((const*ml*cl)/w**2)+(const*mb*sb*sl)/w**2
    
    #V-ml
    jac3[2,1]=(const*cl)/w
    
    #V-mb
    jac3[2,2]=-((const*sb*sl)/w)
    
    #V-vrad
    jac3[2,3]=cb*sl
    
    
    #W-par
    jac3[3,0]=-((const*mb*cb)/w**2)
    #W-ml
    jac3[3,1]=0
    #W-mb
    jac3[3,2]=(const*cb)/w
    #W-vrad
    jac3[3,3]=sb        
    
    return np.dot(jac3,np.dot(jac2,jac0))
    
