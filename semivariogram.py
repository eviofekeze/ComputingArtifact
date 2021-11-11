class semivariogram:
    """

    """
    def solve(self, dist, acc):
        """
        Implements a solver to compute, semi variance, autocorrelation,
        for equally spaced data
        
        variance of a parameter of interest over a space
        
        Input:
            dist = horizontal/vertical distance of value we want to compute for
            acc = value we want to investigate it spatial variability
        Output:
            semi_var = semivariance
            cov_ = variance
            autocorr_ = autocorrelation
        
        """
        
        self.dist = dist
        self.acc = acc
        
        self.dx = np.mean(np.diff(self.dist)) #The average spacing between two points in the space 
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
        
    def plot_variogram(self,dist, acc, figsize = (15,6), color = 'k'):
        temp = semivariogram()
        temp.solve(dist,acc)
        plt.figure(figsize = figsize)
        plt.scatter(temp.h, temp.semi_var_,color = color)
        plt.title('Semivariance Vs lag') # Set Title
        plt.ylabel('Semivariance') # label y axis
        plt.xlabel('Lag (m)')
        
    def plot_covariance(self,dist,acc,figsize = (15,6), color = 'k'):
        temp = semivariogram()
        temp.solve(dist,acc)
        plt.figure(figsize = figsize)
        plt.scatter(temp.h, temp.cov_, color = color) 
        plt.title('Covariance Vs Lag') # Set Title
        plt.ylabel('Covariance') # label y axis
        plt.xlabel('Lag (m)') # # labelx axis
    
    
    def plot_autocorr(self, dist, acc, figsize = (15,6), color = 'k'):
        temp = semivariogram()
        temp.solve(dist,acc)
        plt.figure(figsize = figsize)
        plt.scatter(temp.h, temp.autocorr_, color = color) 
        plt.title('Autocorrelation Vs Lag') # Set Title
        plt.ylabel('Autocorrelation') # label y axis
        plt.xlabel('Lag (m)') # # labelx axis