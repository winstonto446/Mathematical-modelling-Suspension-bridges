import numpy as np
import matplotlib.pyplot as plt
from vibrations_base import VibrationsBase

class SuspensionBridgeBase(VibrationsBase):
    """ Class which computes the normal modes of a suspenssion bridge."""

    def __init__(self, N, theta_tower=30):
        """ Initialiser.
            N : number of nodes
            theta_tower : value of theta at the tower
        """
        self.N=N # number of nodes including the towers
        self.Ksc = 5e7
        self.Kbend = 0
        self.M=15
        
        # angle of each bridge segment from the horizontal 
        self.thetas = np.zeros(self.N, dtype='float64')
        self.set_thetas(theta_tower)

        # Create the vibration matrix. Done in sub classes
        A = self.create_A()
 
        super().__init__(A)

    def create_A(self):
        A = np.identity(self.N)
        return(A)
        
    def set_thetas(self, theta_tower = 30):
       """ Initialise the values of theta
           theta_tower : the angle at the tower
       """
       self.MoT = 2*np.tan(theta_tower*np.pi/180)/(self.N-2)
       midn = self.N//2
       self.thetas[midn] = 0
       for i in range(1, midn):
           self.thetas[midn-i] = -np.arctan(i*self.MoT)
           self.thetas[midn+i] = np.arctan(i*self.MoT)

      
    def plot_spectrum(self, fname=""):   
        """ Plot the vibration spectrum of the bridge in blue.
            Plot multiples of the fundamental frequency in red.
            fname:  Filename for figure output (if non-empty).
        """
        Nmax = self.N-2
        # Plot actual frequency ration in blue
        plt.plot(self.nu, "b*")

        # Plot harmonic frequencies in red
        plt.plot(range(1, Nmax+1)*self.nu[0], "r-")
        plt.xlabel("i", fontsize=24)
        plt.ylabel(r'$\nu$', fontsize=24)
        plt.tight_layout(pad=0.4)
        if(fname != ""):
           plt.savefig(fname)
        plt.show()
