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
    
 
def pmradec2lb(a,d,l,b,mua,mud):
    siphi=np.cos(deo)*np.sin(a-alpm)/np.cos(b)
    cophi=(np.sin(deo)-np.sin(d)*np.sin(b))/(np.cos(d)*np.cos(b))
    
    ml=mua*cophi+mud*siphi
    mb=-mua*siphi+mud*cophi
    
    return ml,mb

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
    l,b=raddec2lb(a,d)

    #trigonometric functions
    
    
    #mas -> rad
    jac0[0,0]=np.pi/180/3600
    jac0[1,1]=np.pi/180/3600
    jac0[2,2]=1
    jac0[3,3]=1
    jac0[4,4]=1
    jac0[5,5]=1
    
    #from equatorial to galactic coordinates
    
    #l-alpha
    jac1[0,0]=((np.cos(incl)/(np.cos(a-alom))**2+np.sin(incl)/np.cos(a-alom)*np.tan(a-alom)*np.tan(d))/(
        1+(np.cos(incl)*np.tan(a-alom)+np.sin(incl)/np.cos(a-alom)*np.tan(d))**2))
    #l-delta
    jac1[0,1]=(np.sin(incl)/(np.cos(a-alom)*(np.cos(d))**2)/(
        1+(np.cos(incl)*np.tan(a-alom)+np.sin(incl)/np.cos(a-alom)*np.tan(d))**2))
    #b-alpha
    jac1[1,0]=(-((np.cos(a-alom)*np.cos(d)*np.sin(incl))/np.sqrt(
        1-(np.cos(incl)*np.sin(d)-np.cos(d)*np.sin(a-alom)*np.sin(incl))**2)))
    
    #b-delta
    jac1[1,1]=((np.cos(d)*np.cos(incl)+np.sin(a-alom)*np.sin(d)*np.sin(incl))/np.sqrt(
        1-(np.cos(incl)*np.sin(d)-np.cos(d)*np.sin(a-alom)*np.sin(incl))**2))
        
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
    jac2[3,0]=-((mud*np.cos(deo)*np.cos(l-theta))/(1-(np.cos(b)*np.cos(deo)*np.cos(
    l-theta)+np.sin(b)*np.sin(
    deo))**2)**(0.5))-(mua*np.cos(deo)*(np.cos(b)*np.cos(deo)*np.cos(
    l-theta)+np.sin(b)*np.sin(deo))*(np.sin(deo)-np.sin(b)*(np.cos(b)*np.cos(
    deo)*np.cos(l-theta)+np.sin(b)*np.sin(deo)))*np.sin(l-theta))/(
    1-(np.cos(b)*np.cos(deo)*np.cos(l-theta)+np.sin(b)*np.sin(deo))**2)**(
    1.5)+(mua*np.cos(deo)*np.sin(b)*np.sin(l-theta))/(
    1-(np.cos(b)*np.cos(deo)*np.cos(l-theta)+np.sin(b)*np.sin(deo))**2)**(
    0.5)+(mud*np.cos(b)*np.cos(deo)**2*(np.cos(b)*np.cos(deo)*np.cos(
    l-theta)+np.sin(b)*np.sin(deo))*np.sin(l-theta)*np.sin(l-theta))/(
    1-(np.cos(b)*np.cos(deo)*np.cos(l-theta)+np.sin(b)*np.sin(deo))**2)**(1.5)
    
    #ml-b
    jac2[3,1]=(mua/(np.cos(b))*(-(np.cos(deo)*np.cos(l-theta)*np.sin(
        b))+np.cos(b)*np.sin(deo))*(np.cos(b)*np.cos(deo)*np.cos(
        l-theta)+np.sin(b)*np.sin(deo))*(np.sin(deo)-np.sin(b)*(np.cos(
        b)*np.cos(deo)*np.cos(l-theta)+np.sin(b)*np.sin(deo))))/(
        1-(np.cos(b)*np.cos(deo)*np.cos(l-theta)+np.sin(b)*np.sin(
        deo))**2)**(1.5)+(mua/(np.cos(b))*(-(np.sin(b)*(
        -(np.cos(deo)*np.cos(l-theta)*np.sin(b))+np.cos(b)*np.sin(
        deo)))-np.cos(b)*(np.cos(b)*np.cos(deo)*np.cos(l-theta)+np.sin(
        b)*np.sin(deo))))/(1-(np.cos(b)*np.cos(deo)*np.cos(l-theta)+np.sin(
        b)*np.sin(deo))**2)**(0.5)-(mud*np.cos(deo)*(-(np.cos(deo)*np.cos(
        l-theta)*np.sin(b))+np.cos(b)*np.sin(deo))*(np.cos(b)*np.cos(
        deo)*np.cos(l-theta)+np.sin(b)*np.sin(deo))*np.sin(l-theta))/(
        1-(np.cos(b)*np.cos(deo)*np.cos(l-theta)+np.sin(b)*np.sin(deo))**2)**(
        1.5)+(mua/(np.cos(b))*(np.sin(deo)-np.sin(b)*(np.cos(b)*np.cos(
        deo)*np.cos(l-theta)+np.sin(b)*np.sin(deo)))*np.tan(b))/(1-(np.cos(
        b)*np.cos(deo)*np.cos(l-theta)+np.sin(b)*np.sin(deo))**2)**(0.5)
    
    #ml-mua
    jac2[3,3]=(1/(np.cos(b))*(np.sin(deo)-np.sin(b)*(np.cos(b)*np.cos(
        deo)*np.cos(l-theta)+np.sin(b)*np.sin(deo))))/(1-(np.cos(b)*np.cos(
        deo)*np.cos(l-theta)+np.sin(b)*np.sin(deo))**2)**(0.5)
    
    #ml-mud
    jac2[3,4]=-((np.cos(deo)*np.sin(l-theta))/(1-(np.cos(b)*np.cos(
        deo)*np.cos(l-theta)+np.sin(b)*np.sin(deo))**2)**(0.5))
    
    
    #mb-l
    jac2[4,0]=(mua*np.cos(deo)*np.cos(l-theta))/(1-(np.cos(b)*np.cos(
        deo)*np.cos(l-theta)+np.sin(b)*np.sin(deo))**2)**(0.5)-(mud*np.cos(
        deo)*(np.cos(b)*np.cos(deo)*np.cos(l-theta)+np.sin(b)*np.sin(
        deo))*(np.sin(deo)-np.sin(b)*(np.cos(b)*np.cos(deo)*np.cos(
        l-theta)+np.sin(b)*np.sin(deo)))*np.sin(l-theta))/(1-(np.cos(b)*np.cos(
        deo)*np.cos(l-theta)+np.sin(b)*np.sin(deo))**2)**(1.5)+(mud*np.cos(
        deo)*np.sin(b)*np.sin(l-theta))/(1-(np.cos(b)*np.cos(deo)*np.cos(
        l-theta)+np.sin(b)*np.sin(deo))**2)**(0.5)-(mua*np.cos(b)*np.cos(
        deo)**2*(np.cos(b)*np.cos(deo)*np.cos(l-theta)+np.sin(b)*np.sin(
        deo))*np.sin(l-theta)*np.sin(l-theta))/(1-(np.cos(b)*np.cos(deo)*np.cos(
        l-theta)+np.sin(b)*np.sin(deo))**2)**(1.5)
    
    #mb-b
    jac2[4,1]=(mud/(np.cos(b))*(-(np.cos(deo)*np.cos(l-theta)*np.sin(b))+np.cos(
        b)*np.sin(deo))*(np.cos(b)*np.cos(deo)*np.cos(l-theta)+np.sin(b)*np.sin(
        deo))*(np.sin(deo)-np.sin(b)*(np.cos(b)*np.cos(deo)*np.cos(
        l-theta)+np.sin(b)*np.sin(deo))))/(1-(np.cos(b)*np.cos(deo)*np.cos(
        l-theta)+np.sin(b)*np.sin(deo))**2)**(1.5)+(mud/(np.cos(b))*(-(
        np.sin(b)*(-(np.cos(deo)*np.cos(l-theta)*np.sin(b))+np.cos(b)*np.sin(
        deo)))-np.cos(b)*(np.cos(b)*np.cos(deo)*np.cos(l-theta)+np.sin(
        b)*np.sin(deo))))/(1-(np.cos(b)*np.cos(deo)*np.cos(l-theta)+np.sin(
        b)*np.sin(deo))**2)**(0.5)+(mua*np.cos(deo)*(-(np.cos(deo)*np.cos(
        l-theta)*np.sin(b))+np.cos(b)*np.sin(deo))*(np.cos(b)*np.cos(
        deo)*np.cos(l-theta)+np.sin(b)*np.sin(deo))*np.sin(l-theta))/(
        1-(np.cos(b)*np.cos(deo)*np.cos(l-theta)+np.sin(b)*np.sin(deo))**2)**(
        1.5)+(mud/(np.cos(b))*(np.sin(deo)-np.sin(b)*(np.cos(b)*np.cos(deo)*np.cos(
        l-theta)+np.sin(b)*np.sin(deo)))*np.tan(b))/(1-(np.cos(b)*np.cos(deo)*np.cos(
        l-theta)+np.sin(b)*np.sin(deo))**2)**(0.5)
    
    #mb-mua
    jac2[4,3]=(np.cos(deo)*np.sin(l-theta))/(1-(np.cos(b)*np.cos(deo)*np.cos(l-theta)+np.sin(b)*np.sin(deo))**2)**(0.5)
         
    #mb-mud
    jac2[4,4]=(1/(np.cos(b))*(np.sin(deo)-np.sin(b)*(np.cos(b)*np.cos(
        deo)*np.cos(l-theta)+np.sin(b)*np.sin(deo))))/(1-(np.cos(b)*np.cos(
        deo)*np.cos(l-theta)+np.sin(b)*np.sin(deo))**2)**(0.5)
    
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
    jac3[3,0]=-((const*ml*np.cos(l))/w)-vrad*np.cos(b)*np.sin(l)+(const*mb*np.sin(b)*np.sin(l))/w
    #U-b
    jac3[3,1]=-((const*mb*np.cos(b)*np.cos(l))/w)-vrad*np.cos(l)*np.sin(b)
    
    #U-par
    jac3[3,2]=(const*mb*np.cos(l)*np.sin(b))/w**2+(const*ml*np.sin(l))/w**2
    
    #U-ml
    jac3[3,3]=-((const*np.sin(l))/w)
    
    #U-mb
    jac3[3,4]=-((const*np.cos(l)*np.sin(b))/w)
    
    #U-vrad
    jac3[3,5]=np.cos(b)*np.cos(l)
    
    
    #V-l
    jac3[4,0]=vrad*np.cos(b)*np.cos(l)-(const*mb*np.cos(l)*np.sin(b))/w-(const*ml*np.sin(l))/w
         
    #V-b
    jac3[4,1]=-((const*mb*np.cos(b)*np.sin(l))/w)-vrad*np.sin(b)*np.sin(l)
    
    #V-par
    jac3[4,2]=-((const*ml*np.cos(l))/w**2)+(const*mb*np.sin(b)*np.sin(l))/w**2
    
    #V-ml
    jac3[4,3]=(const*np.cos(l))/w
    
    #V-mb
    jac3[4,4]=-((const*np.sin(b)*np.sin(l))/w)
    
    #V-vrad
    jac3[4,5]=np.cos(b)*np.sin(l)
    
    
    #W-l
    jac3[5,0]=0
    #W-b
    jac3[5,1]=vrad*np.cos(b)-(const*mb*np.sin(b))/w
    #W-par
    jac3[5,2]=-((const*mb*np.cos(b))/w**2)
    #W-ml
    jac3[5,3]=0
    #W-mb
    jac3[5,4]=(const*np.cos(b))/w
    #W-vrad
    jac3[5,5]=np.sin(b)
    
    
    return np.dot(jac3,np.dot(jac2,np.dot(jac1,jac0)))
