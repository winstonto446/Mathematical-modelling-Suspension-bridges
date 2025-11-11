import numpy as np
import matplotlib.pyplot as plt 
import ode_rk4

class Bridge_base(ode_rk4.ODE_RK4):
   """ A class to define the general attributes of suspended bridges
   """
   def __init__(self, V0=[0], dt=0.1, t0=0, theta_tower=30):
       """ Initialiser
          :param V0 : initial value (as an array or a list) 
          :param dt : integration time step
          :param t0 : initial time
          :param theta_tower: bridge angle at the tower
       """
       super().__init__(V0, dt, t0)
       self.N = len(V0)//4 # Number of nodes
       self.l = 1 # length of each segment
       self.M = 15 # Mass of bridge segment
       self.Ksc =  2e7 # Spring coefficent of suspension cable segment
       self.Kbend = 1e4 # Bending coefficent of suspension cable
       self.thetas = np.zeros(self.N, dtype='float64')
       self.set_thetas(theta_tower)
       self.MoT = 1 #  Ratio M/T for tension in the cables
       self.Gamma = 2 #  Friction coef.
      
   def set_thetas(self, theta_tower = 30):
       """ Initialise the values of theta
           theta_tower : the angle at the tower
           : param theta_tower : bridge angle at the tower
       """
       self.MoT = 2*np.tan(theta_tower*np.pi/180)/(self.N-2)
       midn = self.N//2
       self.thetas[midn] = 0
       for i in range(1, midn):
           self.thetas[midn-i] = -np.arctan(i*self.MoT)
           self.thetas[midn+i] = np.arctan(i*self.MoT)
           
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
      I = np.index_exp[1:self.N-1]
      Im1 = np.index_exp[0:self.N-2]
      Ip1 = np.index_exp[2:self.N]

      L = self.Load(t)
      # IMPROVE SPEED BY USING NO LOOPS 

      for i in range(1,self.N-1):
        FX[i] = VX[i]
        FZ[i] = VZ[i]
        FVX[i] = (self.Ksc*(((Z[i-1]-Z[i])*np.sin(self.thetas[i])\
                           +(X[i-1]-X[i])*np.cos(self.thetas[i]))\
                           *np.cos(self.thetas[i])\
                        +((Z[i+1]-Z[i])*np.sin(self.thetas[i+1])
                         +(X[i+1]-X[i])*np.cos(self.thetas[i+1]))
                           *np.cos(self.thetas[i+1]))- self.Gamma*VX[i])/self.M
        FVZ[i] = (self.Ksc*(((Z[i-1]-Z[i])*np.sin(self.thetas[i])\
                           +(X[i-1]-X[i])*np.cos(self.thetas[i]))\
                           *np.sin(self.thetas[i])\
                        +((Z[i+1]-Z[i])*np.sin(self.thetas[i+1])
                         +(X[i+1]-X[i])*np.cos(self.thetas[i+1]))
                           *np.sin(self.thetas[i+1]))- self.Gamma*VZ[i])/self.M
        FVZ[i] += (self.Kbend*(0.5*(Z[i-1]+Z[i+1])-Z[i]) + L[i])/self.M
        
      return(eq)

   def Load(self, t):
       """ Return a time dependant load which is then added to the force in F
           :param t : current time
       """
       f_load = np.zeros(self.N, dtype='float64')
       # No load
       f_load[self.N//2-1] = f_load[self.N//2] = self.central_load
       return(f_load)

   def plot_XZ_profile(self, t):
       """ Plot the XZ profile of the bridge
           :param t : current time        
       """
       i = 0
       # look for nearest index for the specified time
       while (i < len(self.t_list) and self.t_list[i] < t):
          i += 1
       if (i >= len(self.t_list) ): # Use last index if t too large
          i = len(self.t_list) -1
       X = self.V_list[i][:self.N]
       Z= self.V_list[i][self.N:2*self.N]
       plt.title("t={:.3g}".format(t))
       plt.plot(range(self.N), X, "b*", label="X")
       plt.plot(range(self.N), Z, "r*", label="Z")
       plt.legend(loc='lower right', bbox_to_anchor=(1.20, 0.0), prop={'size':10})
       ax = plt.gca()
       ax.set_aspect(1.0/ax.get_data_ratio()*0.5)
       plt.tight_layout(rect=[0, 0, 1.1, 1])
       plt.xlabel("i", fontsize=24)
       plt.ylabel("X/Z", fontsize=24)
       plt.tight_layout(pad=0.5)

   def plot_trace(self, n, XorZ='Z' , fmt=["k-","r-","g-","b-","c-","m-","y-"]):
      """ Plot the bridge along the X or Z axis ate different times
           :param n : number of profiles
           :param XorZ :  'X' or 'Z' axis.
           :param fmt : plot format to use.
      """
      nt = len(self.t_list)
      step_t = nt//n
      x = np.arange(self.N)
      dx = 0.25*self.N/nt

      shift = 0
      if(XorZ=='Z'): shift = self.N

      # find maximum value
      Vmax = 0
      for i in range(0,nt,step_t):
         m = np.max(np.abs(self.V_list[i][shift:self.N+shift]))
         Vmax = max(m,Vmax)
         

      plt.title(XorZ)
      dv = 2*Vmax/nt
      #print(dx,dv)
      j = 0
      for i in range(0,nt,step_t):
         plt.plot(x+i*dx,self.V_list[i][shift:self.N+shift]+i*dv, fmt[j])
         j += 1
         if(j >= len(fmt)): j = 0
         


         
