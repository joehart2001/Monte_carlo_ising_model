import numpy as np

class IsingLattice:

    E = 0.0
    E2 = 0.0
    M = 0.0
    M2 = 0.0

    n_cycles = 0

    def __init__(self, n_rows, n_cols):
        self.n_rows = n_rows
        self.n_cols = n_cols
        self.lattice = np.random.choice([-1,1], size=(n_rows, n_cols))
        #self.N = n_rows*n_cols
        self.J = 1
        
        # def f(x, a, b, c):
        #     return a + b*np.exp(c*x)
        
        def f(x, a, b):
            return a*np.exp(b*x)
        
        #self.cutoff = f(n_rows, a = -2.00899300e+03, b = 8.71557693e+02, c = 1.59759578e-01) + 15000
        self.cutoff = f(n_rows, a = 10**(-0.52), b=4.11)
    
    def energy(self):
        """
        Return the total energy of the current lattice configuration.
        """
        # shift the lattice by 1 in each direction
        # Just 2 directions to avoid double counting
        R = np.roll(self.lattice, shift=-1, axis=0)
        UP = np.roll(self.lattice, shift=-1, axis=1)
        
        # Compute the sum of neighbouring spins
        sum_sisj = np.sum(self.lattice * (R + UP), axis=(0, 1))
        
        # Compute and return the energy
        energy = -self.J * sum_sisj
        return energy

    def magnetisation(self):
        "Return the total magnetisation of the current lattice configuration."
        magnetisation = np.sum(self.lattice)
        
        return magnetisation

    
        

    def montecarlostep(self, T):
        self.n_cycles += 1
        #can do this in 1-2 if statements
        #generating a random number lets you test the probability distribution for one value
        
        # complete this function so that it performs a single Monte Carlo step
        energy = self.energy()
        
        #the following two lines will select the coordinates of the random spin for you
        random_i = np.random.choice(range(0, self.n_rows))
        random_j = np.random.choice(range(0, self.n_cols))
        
        # flip the spin
        self.lattice[random_i, random_j] = -1*self.lattice[random_i, random_j]

        # calc energy of new config
        new_energy = self.energy()
        E_diff = new_energy - energy
        
        #the following line will choose a random number in the rang e[0,1) for you
        random_number = np.random.random()
        
        botlzmann = np.exp(-E_diff/(T)) #*1.38e10-23
        
        if (E_diff > 0 and random_number > botlzmann):
            # reject
            self.lattice[random_i, random_j] = -1*self.lattice[random_i, random_j]
               
        else:
            # accept the copy and return the new config and magnetisation
            energy = new_energy
        
        M = self.magnetisation()
        
        
        if self.n_cycles > self.cutoff:
            self.E += energy
            self.E2 += energy**2
            self.M += M
            self.M2 += M**2
        else:
            pass
        
        return energy, M

            



    def statistics(self):
        #complete this function so that it calculates the correct values for the averages of E, E*E (E2), M, M*M (M2), and returns them with Nsteps
        #Z = (2*np.cosh(self.J/1))**(self.N)
        N = self.n_cycles
        cutoff = self.cutoff
        if N > cutoff:
            avg_E = self.E/(N-cutoff)
            avg_E2 = self.E2/(N-cutoff)
            avg_M = self.M/(N-cutoff)
            avg_M2 = self.M2/(N-cutoff)
        
        return avg_E, avg_E2, avg_M, avg_M2, N
