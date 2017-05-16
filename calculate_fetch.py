import numpy as np
from scipy.interpolate import RectBivariateSpline, interp1d, interp2d
import soda.dataio.datadownload.getNOAAWeatherStation as noaa
import utm
from init_builder import *
import pickle
import os
genericbathy  = working_dir + '/msl1k.asc'
g = 9.81  # gravity constant
max_fetch_distance = 60000  # maximum fetch distance


class FetchModel:
    def __init__(self, bathyfile=genericbathy, xllcorner=531320.83645285, yllcorner=4142047.4765937, cellsize=100,
                 cellvalue=0.01):
        """
        Initialize the FetchModel class.

        :param bathyfile: Bathymetry (default is "msl1k.asc")
        :param xllcorner: easting of lower left corner
        :param yllcorner: northing of lower left corner
        :param cellsize: horizontal extent of grid scale in meters
        :param cellvalue: depth scaling factor in meters
        """
        x, y, bath = self._readbathy(bathyfile, xllcorner, yllcorner, cellsize, cellvalue)
        self.x = x
        self.y = y
        self.bath = bath
        self.xllcorner = xllcorner
        self.yllcorner = yllcorner
        self.cellsize = cellsize
        self.data = {}
        b =np.rot90(self.bath, 3)

        self.f = RectBivariateSpline(self.x, self.y, b)
    # f = RectBivariateSpline(self.y, self.x, self.bath.transpose())


    @staticmethod
    def _readbathy(bathyfile=genericbathy, xllcorner=531320.83645285, yllcorner=4142047.4765937, cellsize=100,
                   cellvalue=0.01):

        """ Create a grid of bathymetry with extent defined by the bathymetry file and cell size and easting/northing.

        :param bathyfile: Bathymetry (default is "msl1k.asc")
        :param xllcorner: easting of lower left corner
        :param yllcorner: northing of lower left corner
        :param cellsize: horizontal extent of grid scale in meters
        :param cellvalue: depth scaling factor in meters
        :return: x, y, depth
        """
        bath = np.loadtxt(bathyfile) * cellvalue
        bath = bath
        nrows, ncols = bath.shape
        x = xllcorner + cellsize * np.arange(0, ncols, 1)
        y = yllcorner + cellsize * np.arange(0, nrows, 1)

        return x, y, bath

    def _findfetch(self, x0, y0, angle):
        """
        Find the fetch given a wind angle and easting/northing location.

        :param x0: Easting of target point
        :param y0: Northing of target point
        :param angle: Wind angle (angle that the wind is coming from)
        :return: fetch (distance in meters)
        """
        N = 5000
        r = np.linspace(0, max_fetch_distance, N)
        xnew = x0 + np.cos(angle) * r
        ynew = y0 + np.sin(angle) * r
#        print(xnew,ynew)
#        print(self.f(xnew, ynew)[0,:])
#        h = self.f(xnew, ynew)[0,:]
#        # choose the maximum fetch distance
#        if sum(h<=0) == 0:
#            return max_fetch_distance
#        return np.sqrt(np.power(xnew[h<=0][0] - x0, 2.0) + np.power(ynew[h<=0][0] - y0, 2.0))
#        
        for i in range(0,N):
            hq = self.f(xnew[i], ynew[i])[0][0]
#            print(hq)
            if hq <= 0:
                return r[i]
        return max_fetch_distance

    @staticmethod
    def _waveheight(u10, h0, fetch):
        """
        Compute wave height using Army Corp of Engineer's fetch limited wave model

        :param u10: 10 m wind speed
        :param h0: depth of target location
        :param fetch: fetch of target location
        :return: T, Hs (wave period, wave height)
        """
        ua = 0.71 * np.power(u10, 1.23)
        h0[h0 < 0.005] = np.nan*len([h0 < 0.005] )
        froude = g * h0 / ua / ua
        c0 = np.tanh(0.530 * np.power(froude, 3.0 / 4.0))
        c1 = np.tanh(0.833 * np.power(froude, 3.0 / 4.0))
        Hs = np.power(ua, 2) / g * 0.283 * c0 * np.tanh(0.00565 * (np.sqrt(g * fetch / ua / ua)) / c0)
        T = ua / g * 7.54 * c1 * np.tanh(0.0379 * np.power(g * fetch / ua / ua, 1.0 / 3.0) / c1)

        return T, Hs

    def getwaveheight(self, x0, y0, angle, u10, fetch=None):
        """
        Call function for the model that provides the wave period and height at a specified location
        :param x0: Easting of target point
        :param y0: Northing of target point
        :param angle: Wind angle (angle that the wind is coming from aka direction upwind)
        :param u10: 10 m wind speed
        :param fetch: fetch of target location
        :return: T, Hs (wave period, wave height)
        """
        if fetch == None:
            fetch = self._findfetch(x0, y0, angle)

        h0 = self.f(x0, y0)[0,:]
        T, Hs = self._waveheight(u10, h0, fetch)
        return T, Hs, fetch

    @staticmethod
    def _convertwinddata(data):

        # Load in met station data Time in units of minutes since 1970-01-01 00:00:00, lat and long in lat and long, elev in
        # meters, wind speed in m/s.
        """
        Converts the soda weather data to an internally used structure

        :param data: soda data from getNOAAWeatherStation
        :rtype: noaa wind data structure
        """
        time = 0
        data_by_station = {}
        for d in data:
            for k in d.keys():
                if k == 'Uwind' or k == 'Vwind':
                    if not d[k]['StationID'] in data_by_station:
                        x = utm.from_latlon(d[k]['coords'][1]['Value'], d[k]['coords'][0]['Value'])
                        datapoint = {'easting': x[0], 'northing': x[1],
                                     'elev': d[k]['coords'][2]['Value'], 'time': d[k]['coords'][3]['Value']}
                        time = d[k]['coords'][3]['Value'][10]
                    datapoint[k] = d[k]['Data']
                    f = interp1d(datapoint['time'], datapoint[k])
                    datapoint[k+'f'] = f
                    data_by_station[d[k]['StationID']] = datapoint
        return data_by_station

    def downloadwinds(self, timestart = 2016, timeend = 2016):
        """
        Downloads noaa wind data using the soda package getNOAAWeatherStation. This will also convert the data
        structure and save that converted data to self.data

        :rtype: the soda data structure
        """
        ##print(self.xllcorner, self.yllcorner, self.x[-1], self.y[-1])
        llcorner = utm.to_latlon(self.xllcorner, self.yllcorner, 10, 'S')
        urcorner = utm.to_latlon(self.x[-1], self.y[-1], 10, 'S')

        latlon = (llcorner[1], urcorner[1], llcorner[0], urcorner[0])

        localdir = base_dir + 'weather_data/winds'
        ncfile = base_dir + 'weather_data/winds/NCDCNWS_AirObs_'+str(timeend) + '.nc'
        shpfile = base_dir + 'weather_data/winds/NCDCNWS_AirObs_'+str(timeend) + '.shp'
        data = noaa.noaaish2nc(latlon, [timestart, timeend], localdir, ncfile, shpfile)

        self.data = self._convertwinddata(data)
        return data

    def loadwindfromfile(self, filename):
        """
        Loads and converts the soda style and saves as self.data
        :rtype: noaa wind data structure
        """
        with open(filename, 'rb') as inputfile:
            u = pickle._Unpickler(inputfile)
            u.encoding = 'latin1'
            data = u.load()
            # data = pickle.load(inputfile)
        self.data = self._convertwinddata(data)
        return data

    def windinterp(self, x0, y0, time):
        """
        Interpolates
        :param x0: easting  of point of interest
        :param y0: northing of point of interest
        :param time: time of interest
        :return: wind direction (not fetch direction), wind velocity
        """
        x = []
        y = []
        uwind = []
        vwind = []
        for d, k  in zip(self.data.values(),self.data.keys()):
            if k == '720646-99999' or k == '994016-99999':
                continue
            try:
                u = d['Uwindf'](time)
            except:
                u = np.nan
            try:
                v = d['Vwindf'](time)
            except:
                v = np.nan

            if np.isnan(u) or np.isnan(v):
                continue
            uwind.append(u)
            vwind.append(v)
            x.append(d['easting'])
            y.append(d['northing'])
        fuwind = 0
        fvwind = 0

        try:
            fuwind = interp2d(x, y, uwind,kind='linear')
            fvwind = interp2d(x, y, vwind,kind='linear')
        except:
            #print('Error: 2D interpolation Failed, check time and location')
            return np.nan, np.nan
        uwind0 = fuwind(x0, y0)
        vwind0 = fvwind(x0, y0)

        theta = np.arctan2(vwind0, uwind0)
        U10 = np.hypot(uwind0, vwind0)
        return theta, U10, uwind0, vwind0

    # def windinterp(self, x0, y0, time):
    #     """
    #     Interpolates
    #     :param x0: easting  of point of interest
    #     :param y0: northing of point of interest
    #     :param time: time of interest
    #     :return: wind direction (not fetch direction), wind velocity
    #     """
    #     x = []
    #     y = []
    #     uwind = []
    #     vwind = []
    #     for d in self.data.values():
    #         try:
    #             u = d['Uwindf'](time)
    #         except:
    #             u = np.nan
    #         try:
    #             v = d['Vwindf'](time)
    #         except:
    #             v = np.nan
    #
    #         if np.isnan(u) or np.isnan(v):
    #             continue
    #         uwind.append(u)
    #         vwind.append(v)
    #         x.append(d['easting'])
    #         y.append(d['northing'])
    #     fuwind = 0
    #     fvwind = 0
    #
    #     try:
    #         fuwind = interp2d(x, y, uwind)
    #         fvwind = interp2d(x, y, vwind)
    #     except:
    #         #print('Error: 2D interpolation Failed, check time and location')
    #         return np.nan, np.nan
    #     uwind0 = fuwind(x0, y0)
    #     vwind0 = fvwind(x0, y0)
    #
    #     theta = np.arctan2(vwind0, uwind0)
    #     U10 = np.hypot(uwind0, vwind0)
    #     return theta, U10, uwind0, vwind0

    def getwaveheighttime(self, x0, y0, time):
        """
        Returns the waveheight and period for a location and time
        :param x0: point easting
        :param y0: point northing
        :param time: time of interest
        :return: T, Hs (Period and wave height)
        """
        angle, U, u, v = self.windinterp(x0, y0, time)
        #print(angle, U)
        angle = angle[0]-np.pi
        U = U[0]
        T, Hs = self.getwaveheight(x0, y0, angle, U)
        return T, Hs, U, angle, u, v

    @staticmethod
    def savedata(data, filename):
        """
        Saves soda data to a file given as file name
        :param data: soda datastructure
        :param filename: saved data location
        """
        with open(filename, 'wb') as output:
            pickle.dump(data, output, pickle.HIGHEST_PROTOCOL)

    def plotfetch(self, x0, y0, angle, U, time):
        import matplotlib.pyplot as plt
        # plt.imshow(self.bath, extent=[self.x.min(), self.x.max(), self.y.min(), self.y.max()])
        """
        Does some plotting of fetch (mainly for debugging)

        :param x0: easting
        :param y0: northing
        :param angle: fetch direction
        :param U: wind velocity magnitude
        :param time: time of interest (for plotting)
        """
        bathy = self.f(self.x, self.y)
        plt.imshow(np.rot90(bathy), extent=[self.x.min(), self.x.max(), self.y.min(), self.y.max()])
        plt.colorbar()
        plt.scatter(x0, y0, s=50, c='r', marker='o')
        # plt.plot(x0+np.cos(angle-np.pi)*np.linspace(0,1000*U,1000), y0-np.sin(angle-np.pi)*np.linspace(0,10000,1000),c='k')
        plt.plot(x0 + np.cos(angle) * np.linspace(1, 100000, 1000),y0 - np.sin(angle) * np.linspace(1, 100000, 1000),c='k')
        for d in self.data.values():
            plt.scatter(d['easting'], d['northing'], c='m', s=40)
            try:
                u = d['Uwindf'](time)
            except:
                u = np.nan
            try:
                v = d['Vwindf'](time)
            except:
                v = np.nan
            plt.plot(d['easting']+u*np.linspace(0, 1000, 1000), d['northing']+v*np.linspace(0,1000,1000), c='k')
        #print("hi")
        plt.savefig('/Users/jadelson/Desktop/testimage.png')
        #print("here")


# Main function serves as a usage example
if __name__ == "__main__":
    x1 = [5.7000e+05, 5.7002e+05] #easting we care about
    y1 = [4.1577e+06, 4.1597e+06] #northing we care about
    time = 24324168.0 #time we care about

    model = FetchModel() #initalize fetch class
#    data = model.downloadwinds() #download data from noaa website using soda
    #model.savedata(data, '/Users/jadelson/Desktop/data2.pk') #save data to file for saving time later
    # model.loadwindfromfile('/Users/jadelson/Desktop/data2.pk') #load soda data from a file
#    model.savedata(data, '/Users/jadelson/Dropbox/phdResearch/AllOptical/sfbayrsproper/fetch_data.dk')
    # model.savedata(data, '/Users/jadelson/Desktop/data2.pk') #save data to file for saving time later
    fetch_model = FetchModel()
    wind_data_file = base_dir + 'weather_data/weather_data_2014-2017.dk'

    if not os.path.isfile(wind_data_file): 
        wind_station_data = fetch_model.downloadwinds(2014,2017)
        fetch_model.savedata(wind_station_data,wind_data_file)
    else:
        wind_station_data = fetch_model.loadwindfromfile(wind_data_file)
    T, Hs = model._waveheight(np.asarray([10,15]), model.f(x1, y1)[0,:], np.asarray([100, 1000])) #get fetch mode
    
#    #print("start")
#    for x0, y0 in zip(x1, y1):
#        angle, U = model.windinterp(x0, y0, time) #get wind direction and speed
#
#        T, Hs = model.getwaveheighttime(x0, y0, time) #get fetch model wave period and wave height
#    # model.plotfetch(x0, y0, angle+np.pi, U, time)
#
