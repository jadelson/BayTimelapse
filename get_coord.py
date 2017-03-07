import numpy as np
from datetime import datetime
from scipy.optimize import newton


def get_coord(sat_data, lat, lon, m, n, A=None, boxsize=5):
    
    if A is None:
                
        primary = np.array([[sat_data['LON'][0,0],sat_data['LAT'][0,0]],
                    [sat_data['LON'][0,n-1],sat_data['LAT'][0,n-1]],
                    [sat_data['LON'][m-1,n-1],sat_data['LAT'][m-1,n-1]],
                    [sat_data['LON'][m-1,0],sat_data['LAT'][m-1,0]]])
        secondary = np.array([[0, 0],
                      [0, n-1],
                      [m-1, n-1],
                      [m-1., 0]])
        # Pad the data with ones, so that our transformation can do translations too
        n1 = primary.shape[0]
        pad = lambda x: np.hstack([x, np.ones((x.shape[0], 1))])
        unpad = lambda x: x[:,:-1]
        X = pad(primary)
        Y = pad(secondary)
        
    
        # Solve the least squares problem X * A = Y
        # to find our transformation matrix A
        A, res, rank, s = np.linalg.lstsq(X, Y)
    
    transform = lambda x: unpad(np.dot(pad(x), A))
    #A[np.abs(A) < 1e-10] = 0  # set really small values to zero
    newlocal = np.round(np.array([lon, lat, 1]).dot(A))    
    x = int(newlocal[0])
    y = int(newlocal[1])
    #print(newlocal)
    def dist_func(myx,myy,m,n):
        if myx >= m or myx < 0:
#            print("Warning: out of x range")
            return 1000000000.
        if myy >= n or myy < 0:
#            print("Warning: out of y range")
            return 1000000000.
             
        return np.square(sat_data['LAT'][myx,myy]-lat)+np.square(sat_data['LON'][myx,myy]-lon)
    I = 0
    J = 0
    dist = dist_func(x,y,m,n)
    for i in range(-boxsize,boxsize+1):
        for j in range(-boxsize,boxsize+1):
            tempdist = dist_func(x+i,y+j,m,n)
            if tempdist < dist:
                dist = tempdist
                I = i
                J = j
    
    
    #print((dist, x, I, y, J, sat_data['LON'][x+I,y+J], sat_data['LAT'][x+I,y+J]))
    if I > boxsize or J > boxsize:
        print("Warning, boxsize too small!")
    return x+I, y+J, A

