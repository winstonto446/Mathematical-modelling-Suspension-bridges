import numpy as np
import matplotlib.pyplot as plt 
import ode_rk4
from bridge_base import Bridge_base


class Bridge(Bridge_base):
   """ A class to study the general attributes of suspended bridges
       with load attachment
   """
   def __init__(self, V0=[0], dt=0.1, t0=0, theta_tower=30):
       """ Initialiser
          :param V0 : initial value (as an array or a list) 
          :param dt : integration time step
          :param t0 : initial time
          :param theta_tower: bridge angle at the tower
       """
       super().__init__(V0, dt, t0, theta_tower)
       self.central_load = -200
      
   def F(self, t, V):
       """ Evaluate the force on each node at time t.
           : param t : current time
           : param V : displacement of all the nodes as array : [X, Z, VX, VZ,]
       """
       # Array to store the right hand side of the equation
       eq = np.zeros(self.N*4, dtype='float64')
       # A few array refernces to make the code easier to read
       X = V[:self.N]
       Z = V[self.N:2*self.N]
       VX = V[2*self.N:3*self.N]
       VZ = V[3*self.N:]
       
       FX = eq[:self.N]
       FZ = eq[self.N:2*self.N]
       FVX = eq[2*self.N:3*self.N]
       FVZ = eq[3*self.N:]

       cos_th = np.cos(self.thetas)
       sin_th = np.sin(self.thetas)
       

       L = self.Load(t)
       # IMPROVE SPEED BY USING NO LOOPS 

       # perform slicing for the respective elements 
       # to replace the loop: for i in range(1,self.N-1)
       FX[1:self.N-1] = VX[1:self.N-1]
       FZ[1:self.N-1] = VZ[1:self.N-1]
       FVX[1:self.N-1] = (self.Ksc*(((Z[:self.N-2]-Z[1:self.N-1])*sin_th[1:self.N-1]\
                           +(X[:self.N-2]-X[1:self.N-1])*cos_th[1:self.N-1])\
                           *cos_th[1:self.N-1]\
                        +((Z[2:self.N]-Z[1:self.N-1])*sin_th[2:self.N]
                         +(X[2:self.N]-X[1:self.N-1])*cos_th[2:self.N])
                           *cos_th[2:self.N])- self.Gamma*VX[1:self.N-1])/self.M
       FVZ[1:self.N-1] = (self.Ksc*(((Z[:self.N-2]-Z[1:self.N-1])*sin_th[1:self.N-1]\
                           +(X[:self.N-2]-X[1:self.N-1])*cos_th[1:self.N-1])\
                           *sin_th[1:self.N-1]\
                        +((Z[2:self.N]-Z[1:self.N-1])*sin_th[2:self.N]
                         +(X[2:self.N]-X[1:self.N-1])*cos_th[2:self.N])
                           *sin_th[2:self.N])- self.Gamma*VZ[1:self.N-1])/self.M
       FVZ[1:self.N-1] += (self.Kbend*(0.5*(Z[:self.N-2]+Z[2:self.N])-Z[1:self.N-1]) \
                           + L[1:self.N-1])/self.M
        
       return(eq)
   
        

class Bridge_osc(Bridge):
   """ A class to model the oscillation properties of a bridge due to load attachment
   """
   def __init__(self, V0=[0], dt=0.1, t0=0, nu=1, pos=0):
       super().__init__(V0, dt, t0)
       self.nu = nu   # frequency of jumping
       self.pos = pos # index of segment for load
       
   def Load(self, t):
       """ A man jumping on section pos of the bridge at frequency nu
       """
       
       if np.sin(2*np.pi*self.nu*self.t)>0:
           H = np.sin(2*np.pi*self.nu*self.t)
       else:
           H = 0
       # describe the function H(x)
       
       W = -800*H*np.sin(2*np.pi*self.nu*self.t)     
       #the force equation 
       f_load = np.zeros(self.N, dtype='float64')    
       #set the array size
       f_load[self.pos - 1] = f_load[self.pos] = W/4 
       #set the positions of the nodes where the load is applied on
       
       return(f_load)
       
    
   def MaxAmp(self):
       """ Returns the maximum amplitude of vertical displacements generated in 60 seconds
       """
       vertical_displacements = np.abs(np.array(self.V_list)[:, self.N:2*self.N])
       # choose the array with all the rows (all time steps), 
       # and select the respective columns for Z displacements
       # take the absolute value of all the entries selected
       max_amplitude = np.max(vertical_displacements)
       # find the maximum value in the array
       return max_amplitude
       




    
      
       
