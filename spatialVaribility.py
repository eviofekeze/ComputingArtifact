#Helper Functions


import numpy as np
import matplotlib.pyplot as plt
import math


class semivariogram:
    """

    """
    def __init__(self):
        return None
        
    def solve(self, dist, acc):
        """
        Implements a solver to compute, semi variance, autocorrelation,
        for equally spaced data that lie one transect
        
        Input:
            dist = horizontal/vertical distance of value we want to compute for
            acc = value we want to investigate its spatial variability
        Output:
            semi_var = semivariance
            cov_ = variance
            autocorr_ = autocorrelation
        
        """

        self.dist = dist #transect of interest
        self.acc = acc # variable of interest
        self.dx = np.mean(np.diff(self.dist)) #The average spacing between two suceeding points in the space 
        self.extent = np.max(self.dist) - np.min(self.dist) # extent or range of the space
        self.h = np.arange(self.dx,self.extent/2+1,self.dx) # lag distance using half of the extent to avoid bias
        self.N = len(self.dist) # number point data     
        
        self.n_pairs = np.zeros(len(self.dist)) # allocate number of pair     
        self.semi_var_ = np.zeros(len(self.h)) # allocate semi variaince
        self.cov_ = np.zeros(len(self.h)) # allocate variance
        self.autocorr_ = np.zeros(len(self.h)) # allocate autocorrelation 
        
        for i in range(0,len(self.h)):
            self.n_pairs[i] = self.N-i-1 # point data at each lag

            #compute semi variance using formular above
            self.semi_var_[i] =  (1/(2 *self.n_pairs[i])) * np.sum((self.acc[0:(self.N-i-1)] - self.acc[(i+1):self.N])**2)
            
            # compute varinace
            self.cov_[i] = np.mean(self.acc[0:(self.N-i-1)] * self.acc[(i+1):self.N]) - (np.mean(self.acc))**2 #covariance
            
            #compute autocorrelation
            self.autocorr_[i] = self.cov_[i]/np.var(self.acc)
                
        self.semi_var_ = self.semi_var_
        self.cov_ = self.cov_
        self.autocorr_ = self.autocorr_
        
#         if __name__ == '__main__':
#             print('Success!')


class semivariogram2t:
    """

    """
    def __init__(self):
        return None
        
    def solve(self, x_dist, y_dist, var):
        """
        Implements a solver to compute, semi variance for data
        that lie in two transects
        
        Input:
            x_dist = horizontal/vertical distance of value we want to compute for
            y_dist = horizontal/vertical distance of value we want to compute for
            var = value we want to investigate it spatial variability
        Output:
            G = semivariogram
        
        """

        self.x_dist = x_dist # The first transect
        self.y_dist = y_dist # The second transect
        self.var = var # varible of interest
        
        self.dist = np.zeros([len(x_dist), len(self.x_dist)])
        self.gam = np.zeros([len(x_dist), len(self.x_dist)])

        #Compute distance between two pairs of point
        for i in range(0,len(self.x_dist)-1):
            for j in range(i+1, len(self.x_dist)):
                self.dist[i,j] = np.sqrt((self.x_dist[i]-self.x_dist[j])**2 +(self.y_dist[i]-self.y_dist[j])**2) #Euclidean Distance between two points
                self.gam[i,j]= (self.var[i] - self.var[j])**2


        #Fill the other side of the Matrix
        self.dist = self.dist + self.dist.T
        self.gam = self.gam + self.gam.T

        #Lag Distance for estimating the semivariance, we do not want to go beyond half of the maximum distance
        self.hmod = np.arange(0,np.floor(np.max(self.dist)/2),5)

        #Initialse semivariogram vector
        self.G = np.zeros(len(self.hmod))


        #find point within lag range
        for k in range(0,len(self.hmod)-1):
            Ix  = np.logical_and(self.dist >= self.hmod[k], self.dist<= self.hmod[k+1])
            self.G[k] = 0.5*np.mean(self.gam[Ix])
            


class variogramModel:
    """

    """
    def __init__(self):
        return None
    
    def solve(self, h,c,a,type_):
        """
        Variogram using the bounded linear and spherical models
    
        INPUT:
            h = lags
            c = sill
            a = range
            type_ = 'L' for linear and 'S' for spherical
       
        OUTPUT:
            Vm = Modelled varigram
        
        """
        
        self.h = h
        self.c = c
        self.a = a
        self.type_ = type_
        
        
        self.Vm = np.zeros(len(self.h))
        Ix = np.where(self.h < self.a)
        
        
        if (self.type_ == 'L'):
            #lags less than a
            self.Vm[Ix] = self.c*self.h[Ix]/self.a  # Linear 
        elif (self.type_ == 'S'):
            self.Vm[Ix] = self.c*(3*self.h[Ix]/(2*self.a) - 0.5*(self.h[Ix]/self.a)**3) # Spherical 
        
        Ix2 = np.where(self.h>self.a) # lags greater than a
        self.Vm[Ix2] = self.c   
        
        
        
class variogramModelError:
    """

    """
    def __init__(self):
        return None
    
    def solve(self, h,V,c,a,type_):
        """
        Variogram using the bounded linear and spherical models
    
        INPUT:
            h = lags
            V = Experimental Variogram
            c = sill
            a = range
            type_ 'L' for linear and 'S' for spherical
       
        OUTPUT:
            RMSE = The root Mean Squared Error of the modelled values and the Semivariogram
            Vm = The modeled Values
        
        solve(self, h,V,c,a,type_)
        """
        
        self.h = h  # lag distance
        self.V = V # Variogram
        self.c = c # Sill
        self.a = a # Range
        self.type_ = type_ # type S or L
        
        
        self.Vm = np.zeros(len(self.h)) # Preallocate Space for modelled variogram
        Ix = np.where(self.h < self.a) #select index less that range
        
        
        if (self.type_ == 'L'):
            #Solve for index less than Using Linear 
            self.Vm[Ix] = self.c*self.h[Ix]/self.a  # Linear 
        elif (self.type_ == 'S'):
            #Sovel for index less than a using Sperical
            self.Vm[Ix] = self.c*(3*self.h[Ix]/(2*self.a) - 0.5*(self.h[Ix]/self.a)**3) # Spherical 
        
        Ix2 = np.where(self.h>self.a) # lags greater than a
        #Model variogram equals varaince
        self.Vm[Ix2] = self.c
        
        #Compute Root Mean Square
        self.rmse = np.sqrt(np.mean((self.Vm - self.V)**2))
                
        self.Vm = self.Vm
        self.rmse = self.rmse
        
            
def model_variogram(h,c,a,type_):
    """
        Function to compute variogram 
        using either bounded linear or spherical models
    
        INPUT:
            h = lags

            c = sill
            a = range
            type_ = 'L' for linear and 'S' for spherical
       
        OUTPUT:
            V = The modeled Variogram
        
    """
    V = np.zeros(len(h)) # initialize variogram as zeros of lenght h
    Ix = np.where(h <= a) # grab index of all location before the range
    if (type_ == 'L'):
        V[Ix] = c*h[Ix]/a # compute using bounded linear method
    elif(type_ == 'S'):
        V[Ix] = c*(3*h[Ix]/(2*a) - 0.5*(h[Ix]/a)**3) #compute variogram at indexusing bounded spherical method
    
    Ix2 = np.where(h>a) # grab index of locations greater than the range
    V[Ix2] = c # fill greater than the range with calue of the sill
    return V




def model_variogram_error(h,V,c,a,type_):
    """
        Function to computee the square error of a variogram model
    
        INPUT:
            h = lags
            V = Experimental Variogram
            c = sill
            a = range
            type_ = 'L' for linear and 'S' for spherical
       
        OUTPUT:
            RMSE = The root mean squared error of the modelled values of the semivariogram
            V = The modeled Variogram
        
    """
    Vm = np.zeros(len(h)) # initialize variogram as zeros of lenght h
    Ix = np.where(h <= a)  # grab index of all location before the range
    
    if (type_ == 'L'):
        Vm[Ix] = c*h[Ix]/a  #compute variogram at index using bounded linear method
    elif (type_ == 'S'):
        Vm[Ix] = c*(3*h[Ix]/(2*a) - 0.5*(h[Ix]/a)**3) #compute variogram at index using bounded spherical method
        
    Ix2 = np.where(h>a) # grab index of locations greater than the range
    Vm[Ix2] = c # fill greater than the range with calue of the sill
    
    rmse = np.sqrt(np.mean((Vm - V)**2)) #compute RMSE
    
    return rmse


def plot_variogram(dist, acc, figsize = (15,6), color = 'k', label = None):
    """
    Function to plot semivariogram
    
    INPUT:
        dist: Lag distances
        acc: Measurement at lag distance
    OUTPUT:
        Scatter plot of semivariogram 
    """
    temp = semivariogram() # instantiate semivariogram class
    temp.solve(dist,acc) # solve semivariogram
    plt.figure(figsize = figsize) 
    plt.scatter(temp.h, temp.semi_var_,color = color, label = label)
    plt.title('Semivariance Vs lag',fontsize = 16) # Set Title
    plt.ylabel('Semivariance',fontsize = 16) # label y axis
    plt.xlabel('Lag (m)',fontsize = 16)
        
def plot_covariance(dist,acc,figsize = (15,6), color = 'k', label = None):
    """
    Function to plot covariance
    
    INPUT:
        dist: Lag distances
        acc: Measurement at lag distance
    OUTPUT:
        Scatter plot of covariance
    """
    temp = semivariogram() #instantiate semivariogram class
    temp.solve(dist,acc) # solve semivariogram
    plt.figure(figsize = figsize)
    plt.scatter(temp.h, temp.cov_, color = color, label = label) 
    plt.title('Covariance Vs Lag',fontsize = 16) # Set Title
    plt.ylabel('Covariance',fontsize = 16) # label y axis
    plt.xlabel('Lag (m)',fontsize = 16) # # labelx axis
    
    
def plot_autocorr(dist, acc, figsize = (15,6), color = 'k', label = None):
    """
    Function to plot autocorrolation
    
    INPUT:
        dist: Lag distances
        acc: Measurement at lag distance
    OUTPUT:
        Scatter plot of autocorrelation
    """
    temp = semivariogram() #instantiate semivariogram class
    temp.solve(dist,acc) # solve semivariogram
    plt.figure(figsize = figsize)
    plt.scatter(temp.h, temp.autocorr_, color = color, label = label) 
    plt.title('Autocorrelation Vs Lag',fontsize = 16) # Set Title
    plt.ylabel('Autocorrelation',fontsize = 16) # label y axis
    plt.xlabel('Lag (m)',fontsize = 16) # # labelx axis