import numpy as np
import matplotlib.pyplot as plt
from vibrations_base import VibrationsBase

class SuspensionBridge_x(VibrationsBase):

    def create_A(self):
        # TO DO
        pass

class SuspensionBridge_z(SuspensionBridge_x):

     def create_A(self):
        # TO DO
        pass
       
class SuspensionBridge_xz(SuspensionBridge_x):
    def create_A(self):
        self.Kbend = 1e4
        k = self.Kbend/self.Ksc
        # MATRIX FOR XZ BRIDGE
        # IMPROVE BY REMOVING 1 LOOP
        n =  self.N-2
        A = np.zeros([2*n,2*n], dtype='float64')
        # X[i] = V[i]; Z[i] = V[n+i]

        for i in range (1,N-1): # we do not compute A for end nodes
          for j in range (1,N-1):
            if (i==j):
              # Array A must start with index 0 -> shift index by 1
              A[i-1,j-1] = -(np.cos(self.thetas[i])**2+np.cos(self.thetas[i+1])**2)
              # n starts at correct place
              A[n+i-1,n+j-1] = -(np.sin(self.thetas[i])**2+np.sin(self.thetas[i+1])**2) -k
              A[i-1,n+j-1] = -np.cos(self.thetas[i])*np.sin(self.thetas[i])\
                            -np.cos(self.thetas[i+1])*np.sin(self.thetas[i+1])
              A[n+i-1,j-1] = -np.cos(self.thetas[i])*np.sin(self.thetas[i])\
                      -np.cos(self.thetas[i+1])*np.sin(self.thetas[i+1])
            if(i+1 == j):
              A[i-1,j-1] = np.cos(self.thetas[i+1])**2
              A[n+i-1,n+j-1] = np.sin(self.thetas[i+1])**2 +0.5*k 
              A[i-1,n+j-1] = np.cos(self.thetas[i+1])*np.sin(self.thetas[i+1])
              A[n+i-1,j-1] = np.cos(self.thetas[i+1])*np.sin(self.thetas[i+1])
            if (i-1 == j):  
              A[i-1,j-1] = np.cos(self.thetas[i])**2 
              A[n+i-1,n+j-1] = np.sin(self.thetas[i])**2  +0.5*k  
              A[i-1,n+j-1] = np.cos(self.thetas[i])*np.sin(self.thetas[i])
              A[n+i-1,j-1] = np.cos(self.thetas[i])*np.sin(self.thetas[i])
        
        A = A*self.Ksc/self.M
        VibrationsBase.__init__(self,A)
        
