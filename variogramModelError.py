class variogramModelError:
    """

    """
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
            #Solve for place less than Using Linear 
            self.Vm[Ix] = self.c*self.h[Ix]/self.a  # Linear 
        elif (self.type_ == 'S'):
            #Sovel for place less than a using Sperical
            self.Vm[Ix] = self.c*(3*self.h[Ix]/(2*self.a) - 0.5*(self.h[Ix]/self.a)**3) # Spherical 
        
        Ix2 = np.where(self.h>self.a) # lags greater than a
        #Model variogram equals varaince
        self.Vm[Ix2] = self.c
        
        #Compute Root Mean Square
        self.rmse = np.sqrt(np.mean((self.Vm - self.V)**2))
                
        self.Vm = self.Vm
        self.rmse = self.rmse