import numpy as np
import matplotlib.pyplot as plt
from vibrations_base import VibrationsBase
from suspensionbridge_base import SuspensionBridgeBase

class SuspensionBridge_x(SuspensionBridgeBase):
    """ A class which computes the approximated normal modes of 
        the horizontal displacements of a suspension bridge
    """
    def create_A(self):
        """Creates array A corresponding to the equation for
           horizontal displacements of the bridge
        """
        n =  self.N-2
        A = np.zeros([n,n], dtype='float64')                              

        for i in range(1,n):                                              
        # create a loop to set the corresponding values in matrix A
            # if (j==i), then delta_{i,j} = 1
            A[i-1,i-1] = -self.Ksc/self.M*(np.cos(self.thetas[i])**2 \
                          +np.cos(self.thetas[i+1])**2)                   
            # if (j==i+1), then delta_{i,j} = 1
            A[i-1,i] = self.Ksc/self.M*np.cos(self.thetas[i+1])**2      
            # if (j==i-1), then delta_{i,j} = 1
            A[i,i-1] = self.Ksc/self.M*np.cos(self.thetas[i+1])**2        
        # set the entry in A for the last term where i == j == N-2
        A[n-1,n-1] = -self.Ksc/self.M*(np.cos(self.thetas[n])**2 \
                        + np.cos(self.thetas[n+1])**2)                    
        return(A)
        

class SuspensionBridge_z(SuspensionBridge_x):
     """ A class which computes the approximated normal modes of 
         the vertical displacements of a suspension bridge
     """
     def create_A(self):
        """Creates array A corresponding to the equation for
           vertical displacements of the bridge
        """
        n =  self.N-2
        A = np.zeros([n,n], dtype='float64')                     
        for i in range(1,n):                                           
        # create a loop to set the corresponding values in matrix A
            # if (j==i), then delta_{i,j} = 1
            A[i-1,i-1] = -self.Ksc/self.M*(np.sin(self.thetas[i])**2 \
                            + np.sin(self.thetas[i+1])**2)              
            # if (j==i+1), then delta_{i,j} = 1
            A[i-1,i] = self.Ksc/self.M*np.sin(self.thetas[i+1])**2    
            # if (j==i-1), then delta_{i,j} = 1 
            A[i,i-1] = self.Ksc/self.M*np.sin(self.thetas[i+1])**2      
        # set the entry in A for the term where i == j == N-2
        A[n-1,n-1] = -self.Ksc/self.M*(np.sin(self.thetas[n])**2 \
                        + np.sin(self.thetas[n+1])**2)                  
        
        return(A)
       
class SuspensionBridge_xz(SuspensionBridge_x):
    """ A class which computes the normal modes of 
        both vertical and horizontal displacements of a suspension bridge
    """
    def create_A(self):
        """Creates array A corresponding to the equations for both
           horizontal and vertical displacements of the bridge
        """
        self.Kbend = 1e4
        k = self.Kbend/self.Ksc
        # MATRIX FOR XZ BRIDGE
        # IMPROVE BY REMOVING 1 LOOP
        n =  self.N-2
        A = np.zeros([2*n,2*n], dtype='float64')
        # X[i] = V[i]; Z[i] = V[n+i]

        for i in range(1,n): # we do not compute A for end nodes
          # Array A must start with index 0 -> shift index by 1
          # if (i==j), then delta_{i,j}=1
            A[i-1,i-1] = -(np.cos(self.thetas[i])**2+np.cos(self.thetas[i+1])**2)  
            
            A[n+i-1,n+i-1] = -(np.sin(self.thetas[i])**2+np.sin(self.thetas[i+1])**2) -k
            A[i-1,n+i-1] = -np.cos(self.thetas[i])*np.sin(self.thetas[i])\
                          -np.cos(self.thetas[i+1])*np.sin(self.thetas[i+1])
            A[n+i-1,i-1] = -np.cos(self.thetas[i])*np.sin(self.thetas[i])\
                      -np.cos(self.thetas[i+1])*np.sin(self.thetas[i+1])
        
          # if (i+1 == j), then delta_{i,j}=1
            A[i-1,i] = np.cos(self.thetas[i+1])**2                                 
            A[n+i-1,n+i] = np.sin(self.thetas[i+1])**2 +0.5*k 
            A[i-1,n+i] = np.cos(self.thetas[i+1])*np.sin(self.thetas[i+1])
            A[n+i-1,i] = np.cos(self.thetas[i+1])*np.sin(self.thetas[i+1])
            
          # if (i-1 == j), then delta_{i,j}=1 
            A[i,i-1] = np.cos(self.thetas[i+1])**2                                 
            A[n+i,n+i-1] = np.sin(self.thetas[i+1])**2  +0.5*k  
            A[i,n+i-1] = np.cos(self.thetas[i+1])*np.sin(self.thetas[i+1])
            A[n+i,i-1] = np.cos(self.thetas[i+1])*np.sin(self.thetas[i+1])
        
        # set the entries in A for the term where i == j == N-2
        A[n-1,n-1] = -(np.cos(self.thetas[n])**2+np.cos(self.thetas[n+1])**2)      
        # n starts at correct place                                                                     
        A[n+n-1,n+n-1] = -(np.sin(self.thetas[n])**2+np.sin(self.thetas[n+1])**2) -k
        A[n-1,n+n-1] = -np.cos(self.thetas[n])*np.sin(self.thetas[n])\
                      -np.cos(self.thetas[n+1])*np.sin(self.thetas[n+1])
        A[n+n-1,n-1] = -np.cos(self.thetas[n])*np.sin(self.thetas[n])\
                  -np.cos(self.thetas[n+1])*np.sin(self.thetas[n+1])
                  
        A = A*self.Ksc/self.M
        VibrationsBase.__init__(self,A)
        
        return(A)

