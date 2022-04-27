import numpy as np
import matplotlib.pyplot as plt
import random
import math

class Ice(object):
    def __init__(self, N, J=1.0, beta=1, h=0.0, eps=1, T=1):
        """
        : param N : lattice size is N x N
        : param J : interaction strength
        : param h : external magnetic field
        : param T : temperature
        """
        self.N = N
        self.J = J            
        self.h = h
        self.beta = beta
        self.eps = eps
        self.T = T
        self.setup_lattice()
        self.setup_plotting()
        self.Eav = 0
        #self.f = 0
        self.live_show = True

    def setup_lattice(self):
        """ Setup the lattice as a form of a checkerboard. 
        """
        self.lattice = np.ones([self.N, self.N])
        for x in range(0,self.N):
            for y in range(0,self.N):
                self.lattice[x,y]=(x+y)%2
        
    def setup_plotting(self):
        """ Setup a figure for interactive plotting. 
        """
        plt.close("all")
        self.live = plt.figure(figsize=(4, 4))
        plt.ion()
        plt.axis('off')
        plt.show()

    def plot_lattice(self, thermalising=False):
        """ Plot the spin configuration. 
        : param thermalising : if true display lattice in grey instead of blue
        """
        X, Y = np.meshgrid(range(self.N), range(self.N))
        cm = plt.cm.Blues
        if thermalising:
            cm = plt.cm.Greys
        #plt.cla()
        plt.pcolormesh(X, Y, self.lattice, cmap=cm)
        plt.axis('off')
        plt.draw()
        plt.pause(0.01)
        
    
    def monte_carlo(self):
        """Single step of Monte-Carlo using the Cluster-type algorithm
        : inCluster : builds a cluster - zero if the plaquette is in the cluster, and 1 if it is not
        : clusterList : list of points in the cluster
        """
        x = np.random.randint(0, self.N)
        y = np.random.randint(0, self.N)
            
        inCluster = np.zeros([self.N,self.N])
        clusterList = [];
        stack = []; #keep track of which spins to consider (i.e. the nearest neighbours)
        
        Padd = 1 - math.exp(-self.eps/self.T)
            
        colorA = self.lattice[x,y];
        colorB = self.otherColors(colorA)[random.randint(0,1)]
        stack.append((x,y));
            
        #Loop through the stack and consider all the neighbours
        while (len(stack) > 0):
            #print("Stack size: ", len(stack), "x ", x, " y ",y);
            #pick the first item from the stack
            x = stack[0][0]
            y = stack[0][1]
                
            #now that we are considering it, remove it from the stack
            del stack[0]
            
            #ask which color we are currently considering. 
            colorConsidered = self.lattice[x,y]
            if (colorConsidered == colorA or colorConsidered == colorB):
                
                #add it to our global cluster
                inCluster[x,y]=1
                clusterList.append((x,y))
                
                #its neighbours: 
                neighbours = self.returnNeighbours(x,y)
                for q in neighbours:
                    if (inCluster[q[0],q[1]]==0 and not ((q[0],q[1]) in stack)): # i.e. not already in the cluster, and not in the list to consider
                        stack.append((q[0],q[1]))
                
                #add next-nearest neighbours of the same colour with probability Padd < 1
                nextneighbours = self.returnNextNeighbours(x,y)
                if np.random.rand() < Padd:
                    for q in nextneighbours:
                        if colorConsidered == self.lattice[q[0],q[1]]:
                            clusterList.append((q[0],q[1]))
        
        #switch the colors A and B
        for (x,y) in clusterList:
            if (self.lattice[x,y]==colorA):
                self.lattice[x,y]=colorB
            elif (self.lattice[x,y]==colorB):
                self.lattice[x,y]=colorA
            else:
                print('Error!') #this should never happen. 
                        
        sum1 = np.sum(s.lattice==np.roll(s.lattice,2,0))+np.sum(s.lattice==np.roll(s.lattice,-2,0))
        sum2 = np.sum(s.lattice==np.roll(s.lattice,2,1))+np.sum(s.lattice==np.roll(s.lattice,-2,1)) 
        
        E = (self.N*self.eps - self.eps*(sum1+sum2))
        #self.ev = np.var(Es)
        
        return inCluster, clusterList, E
        
    def returnNeighbours(self,x,y):
        """ Return the nearest-neighbours of a given point
        """
        return ((((x+1)%self.N,y),((x-1)%self.N,y),(x,(y+1)%self.N),(x,(y-1)%self.N))); 
    
    def returnNextNeighbours(self,x,y):
        """ Return the next nearest-neighbours of a given point
        """
        return ((((x+2)%self.N,y),((x-2)%self.N,y),(x,(y+2)%self.N),(x,(y-2)%self.N)));

    def otherColors(self,color):
        # this returns the other colors
        if (color == 0):
            return (1,2)
        elif (color == 1):
            return (0,2)
        else:
            return (0,1)
        
    def simulate ( self , temperatures , steps ) :
        """ Do a simulation for a range of values of the temperature .
        : param temperatures : list of temperature to evaluate
        : param steps : number of steps for each Monte Carlo .
        """
        
        # first, a few thermalizing sweeps for smoother data
        for i in range(200*self.N*self.N):
            self.monte_carlo();
            
        results =[]
        Es = []
        E2 = []
        C2 = []
        #self . live_show = False
        l = self . setup_lattice () # init from random
        for T in temperatures :
            self . T = T
            for q in range(1,1000):
                cluster, clusterList, E = self.monte_carlo()
                Es.append(E)
            Eav = np.mean(Es)
            nt = (self.N**2)*(self.T**2)
            C2.append(np.var(Es))
            Cav = np.mean(C2)/nt
            results . append ( (T , Eav, Cav ) )
            print ('T =%12f <E>=%12f <C>=%12f'%( T , Eav, Cav ) )
        return results
        
    def plot_energy(self, T_E_values):
        """
        Plot the magnetisation as a function of the temperature.
        : param T_E_values : list of (T, E)
        """
        plt.close("all")
        plt.ioff()
        plt.axis('on')
        plt.xlabel('Temperature')
        plt.ylabel('Energy')
        plt.title('Energy for eps = 1 and N = 15')
        plt.axvline(x=1/math.log(2), color='m')
        T , E, C = zip (* T_E_values )
        plt . plot (T , E, color='r')
        plt . show ()
    
    def plot_specific(self, T_E_values):
        """
        Plot the magnetisation as a function of the temperature.
        : param T_E_values : list of (T, E)
        """
        plt.close("all")
        plt.ioff()
        plt.axis('on')
        plt.xlabel('Temperature')
        plt.ylabel('Specific Heat Capacity')
        plt.title('Specific Heat Capacity for eps = 0.5 and N = 15')
        T , E, C = zip (* T_E_values )
        plt . plot (T , C, color='k')
        plt . show ()


if __name__ == "__main__":
    
    """
    Simulate and plot energy and specific heat capacity for just N = 20
    """
    
    s = Ice(15)
    s.plot_lattice()
    s.returnNeighbours(0,0)
    s.otherColors(1)
    s. live_show = False
    
    
    print(s.lattice)
   
    results = s . simulate ( np . linspace (0.0 , 4.0 , 30) , 10000 )
    s.plot_lattice()
        
    #s.plot_specific (results)
    s.plot_energy(results)
    
    """
    Plot energies and specific heat capcities for varied N as functions of temperatures
    # """
    # Ns = [10, 15, 20]
    # numNs = len(Ns)
    # numtemperatures = 20
    # temperatures = np.linspace(0.0,4.0,numtemperatures)
    # Eav = np.zeros(shape=(numNs,numtemperatures))
    # Cav = np.zeros(shape=(numNs,numtemperatures))

    # for j,N in enumerate(Ns):
    #     print("Currently on N = ", N)
    #     s = Ice(N)
    #     s.returnNeighbours(0,0)
    #     s.otherColors(1)
    #     s. live_show = False;
    #     Es =[]
    #     E2 =[]
    #     C2 =[]
    #     for i,temp in enumerate(temperatures):
    #         s.T = temp
    #         #print (temp)
    #         for q in range(1,2000):
    #             cluster, clusterList, E = s.monte_carlo()
    #             Es.append(E)
    #             E2.append(E**2)
    #         Eav[j,i] = np.mean(Es)
    #         nt = (N**2)*(s.T**2)
    #         C2.append (np.var(Es))
    #         #print (C2)
    #         Cav[j,i] = np.mean(C2)/nt
    #         #print (Eav)
    #         print (Cav)
    # for j,N in enumerate(Ns):
    #     print (temperatures)
    #     plt.plot(temperatures, Cav[j,:], label=N)
    #     #plt.plot(temperatures, Eav[j,:], label=N)
    #     plt.legend()
    
    # plt.title('Specific Heat Capacity for for N = 10,15,20')
    # #plt.title('Energy for for N = 20,30,40')
    # plt.xlabel('Temperautre')
    # plt.ylabel('Specific Heat Capacity')
    # #plt.ylabel('Energy')
    
