import numpy as np

def radec2lb(ra,dec):
    """
    Returns galactic longitude and latitude (radians) given
    the right ascention and declination (radians)
    """
    #constants
    incl=np.deg2rad(62.87124882)
    alom=np.deg2rad(282.8594813)
    lom=np.deg2rad(32.93680516)
    
    #first part
    sinb=np.sin(dec)*np.cos(incl)-np.cos(dec)*np.sin(incl)*np.sin(ra-alom)
    b=np.arcsin(sinb)
    
    #second part
    sinl=(np.sin(dec)*np.sin(incl)+np.cos(dec)*np.cos(incl)*np.sin(ra-alom))/np.cos(dec)
    cosl=np.cos(dec)*np.cos(ra-alom)/np.cos(dec)
    l=np.arctan2(sinl,cosl)+lom
    if l<0:l=l+2*np.pi
    
    return l,b
    
 
