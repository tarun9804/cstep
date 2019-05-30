#
#distance between two coordinate (lat, lon)
#d = 2*r*asin(sqrt(sin((l2-l1)/2)^2+cos(l2)*cos(l1)*sin((L2-L1)/2)^2))
#r = 6378137m earth radius at equator



import numpy as np
R = 6378137     #radius of earth in meters

#Haversine formula to calculate distance between two coordinate (lat, lon)
def arclength(lat1,lat2,lon1,lon2):
    la1=np.deg2rad(lat1)
    la2=np.deg2rad(lat2)
    lo1=np.deg2rad(lon1)
    lo2=np.deg2rad(lon2)
    h1=np.square(np.sin((la2-la1)/2))
    h2=np.square(np.sin((lo2-lo1)/2))
    h=np.sqrt(h1+np.cos(la1)*np.cos(la2)*h2)
    return 2*R*np.arcsin(h)

#function to find the nearest lat lon from list of pre-defined latitude longitude
def findNearestLatLon(X,Y):
    #find the index of lat lon which is just less than the given lat lons
    i=((lat_array - X)>0).argmax()-1
    j=((lon_array - Y)>0).argmax()-1
    d=np.zeros((4))
    #find the distance from given lat,lon to four corner lat lons
    d[0]=arclength(lat_array[i],X,lon_array[j],Y)
    d[1]=arclength(lat_array[i],X,lon_array[j+1],Y)
    d[2]=arclength(lat_array[i+1],X,lon_array[j],Y)
    d[3]=arclength(lat_array[i+1],X,lon_array[j+1],Y)
    #find the least distance from the four corners
    res=d.argmin()
    if res==0:
        return (lat_array[i],lon_array[j])
    elif res==1:
        return (lat_array[i],lon_array[j+1])
    elif res==2:
        return (lat_array[i+1],lon_array[j])
    else:
        return (lat_array[i+1],lon_array[j+1])



#function to find the nearest lat lon from list of pre-defined latitude longitude
def findLatLon(X,Y):
    #find the index of lat from pre defined list of latitudes which has least difference
    i=(np.abs(lat_array - X)).argmin()
    #find the index of lon from pre defined list of longitudes which has least difference
    j=(np.abs(lon_array - Y)).argmin()
    #return the lat lon with the corresponding index
    return (lat_array[i],lon_array[j])

X=10.76222      #latitude
Y=80.15197     #longitude
lat_array = np.arange(5.05, 37.95, 0.1)  #latitude
lon_array = np.arange(67.05, 97.95, 0.1) #longitude
X = np.linspace(15.05, 17.05, 1000)
Y = np.linspace(77.05, 79.95, 1000)
t=0
#findLatLon(100,100)
for i in range(1000):
    #X=np.random.uniform(5.05, 37.95)
    #Y=np.random.uniform(67.05, 97.95)
    #print(findNearestLatLon(X[i],Y[i]))
    #print(findLatLon(X,Y))
    #print("\n")
    a,b=findNearestLatLon(X[i],Y[i])
    c,d=findLatLon(X[i],Y[i])
    t=t+(a-c)
    t=t+(b-d)
    if t>0:
        print(X[i],Y[i])
        break
    #print(a,b,c,d)
