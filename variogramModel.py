class variogramModel:
    """

    """
    def solve(self, h,c,a,type_):
        """
        Variogram using the bounded linear and spherical models
    
        INPUT:
            h = lags
            c = sill
            a = range
            typre_ 'L' for linear and 'S' for spherical
       
        OUTPUT:
        Vm: Modelled varigram
        
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
        
        
        self.Vm = self.Vm