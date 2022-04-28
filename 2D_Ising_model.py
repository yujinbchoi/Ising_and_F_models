#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 17:44:13 2022

@author: choiyujin
"""

# Here we do a quick simulation of our own 2d ising code.

import numpy as np
import matplotlib.pyplot as plt
import random

class Spins(object):
    def __init__(self, N, J=1, T=1):
        """
        : param N : lattice size is N x N
        : param J : interaction strength
        : param h : external magnetic field
        : param T : temperature
        """
        self.N = N
        self.J = J            
        #self.h = h
        self.T = T
        self.setup_lattice()
        #self.setup_plotting()
        self.m = 0
        self.m2 = 0
        self.m4 = 0


    def setup_lattice(self):
        """ Setup the lattice with a configuration of all-up spins. 
        """
        self.lattice = np.ones([self.N, self.N])
            
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
        #plt.cla()

        
    def mcStep(self):
        """Single step of Monte-Carlo 
        : neighbours : nearest-neighbours of a randomly chosen site on a lattice
        : Ups : the acceptance ratio for the Metropolis algorithm
        """        
        x = np.random.randint(0, self.N)
        y = np.random.randint(0, self.N)
        
        neighbours = self.lattice[(x+1)%self.N,y]+self.lattice[(x-1)%self.N,y]+self.lattice[x,(y+1)%self.N]+self.lattice[x,(y-1)%self.N];
        deltaE = 2*self.lattice[x,y]*neighbours
        Ups = np.exp(-deltaE/self.T)
        
        # flip the spin according to metropolis: 
        if (np.random.rand()<Ups):
            #print("Flipping the spin")
            self.lattice[x,y]=-self.lattice[x,y]
    
    def calcMagnetization(self):
        """ Calculate and return the average magnetisation
        """
        return (np.sum(self.lattice)/(self.N*self.N))
    
    def calcEnergy(self):
        """ Calculate and return the average internal energy
        """        
        sum2 = np.sum(self.lattice*(np.roll(self.lattice, 1, 0) + np.roll(self.lattice, 1, 1)))
        return -self.J*sum2
    
    def mcRun(self,numSweeps,sampleRate=200):
        """ Perform a desired number of sweeps of the lattice 
        and store the average magnetizations in the data structure
        """
        # first, a few thermalizing sweeps for smoother data
        for i in range(200*self.N*self.N):
            self.mcStep();
            
        # initialize the variables
        
        Ms = np.zeros(int((numSweeps*self.N*self.N)/sampleRate));
        Es = np.zeros(int((numSweeps*self.N*self.N)/sampleRate));


        # just do a bunch of steps
        for i in range(numSweeps*self.N*self.N):
            #self.mcStep()
            self.mcStep()
            if (i%(sampleRate) == 0):
                # this variable out just tells us which number of output we are on. 
                out = int(i/(sampleRate));
                #print(out)
                Ms[out]=self.calcMagnetization()
                Es[out]=self.calcEnergy()
                #self.plot_lattice()
                
        # plot the lattice at the end. 
        #self.plot_lattice()
        #plt.plot(Ms[20:])
        
        # store the moments (not the cumulants)
        self.m = np.mean(Ms)
        self.m2 = np.mean(np.power(Ms,2))
        #self.m4 = np.mean(np.power(Ms,4))
        self.mv = np.var(abs(Ms))
        self.ev = np.var(Es)
        self.e = np.mean(Es)
        self.e2 = np.mean(np.power(Es,2))
        
        return Ms


# now initialize and run the model 

# these are the system sizes 
Ns = [40, 60, 80, 100]
numNs = len(Ns)

# these are the values of temperature
numTemperatures = 30
temperatures = np.linspace(1.0,4.0,numTemperatures)

Ms = np.zeros(shape=(numNs,numTemperatures))
sus = np.zeros(shape=(numNs,numTemperatures))
#specific = np.zeros(shape=(numNs,numBetas))

for j,N in enumerate(Ns):
    print("Currently on N = ", N)
    s = Spins(N);
    # now loop through some temperatures
    for i,beta in enumerate(temperatures):
        # set this to a value
        s.T = beta
        s.mcRun(2000)
        Ms[j,i] = abs(s.m)
        #print (Ms)
        # note the M2s is the second cumulant, not the second moment
        #sus[j,i] = (N*s.mv)/s.T
        #sus[j,i] = N*(s.m2 - np.power(s.m,2))/s.T
        print (sus)
        #specific[j,i] = s.ev/(N*s.T*s.T)
        #specific[j,i] = (s.e2 - np.power(s.e,2))/(s.T*s.T*N*N)
        #print (specific)

# now make the plot
for j,N in enumerate(Ns):
    print (temperatures)
    plt.plot(temperatures,Ms[j,:],label=N)
    #plt.legend()

plt.title('Magnetisation for for N = 40,60,80,100')
#plt.title('Energy for for N = 40,60,80,100')
#plt.title('Specific Heat Capacity for for N = 40,60,80,100')
#plt.title('Susceptibility for for N = 40,60,80,100')
plt.xlabel('Temperature')
plt.ylabel('Magnetisation')
#plt.ylabel('Energy')
#plt.ylabel('Specific Heat Capacity')
#plt.ylabel('Susceptibility')

