#!/usr/bin/env python3
__description__ = \
"""
Class for doing a non-linear regression on phage display output data.
"""
__author__ = "Michael J. Harms"
__date__ = "2015-01-20"

import sys, time
import pickle, scipy
import numpy as np

from scipy.optimize import minimize, curve_fit

"""

A = ln(theta_x_0)
B = ln(E_x)

bounds on A: < 0
bounds on B: None
hard part is that sum(theta_0) <= 1.

theta_x(i) = theta_x_0*E_x**i/sum(Q)
           = theta_x_0*exp(i*B)/sum(Q)
           = exp(A)*E_x**i/sum(Q)
           = exp(A + B*i)/sum(Q)

ln(theta_x(i)) = A + B*i - ln(sum(Q))
               = ln(theta_x_0) + B*i - ln(sum(Q))
               = A + ln(E_x)*i - ln(sum(Q))
               = ln(theta_x_0) + ln(E_x)*i - ln(sum(Q))
"""



class fitModel:
    """
    """

    def __init__(self,patterns,degeneracy,rounds_start_at=0):
        """
        Using the data in patterns and degeneracy, create an observable set, 
        objective function, constraints, and bounds taht can then be minimized.
        """

        self.patterns = patterns
        self.degeneracy = degeneracy 
        self.rounds_start_at = rounds_start_at
        
        self.num_patterns = len(patterns)
        self.num_rounds = len(patterns[0,:])
        self.round_exponents = np.array([i + self.rounds_start_at
                                         for i in range(self.num_rounds)])
        self.log_degeneracy = np.log(self.degeneracy)

        # Create observable matrix (thetas x rounds)
        self.y_obs = patterns[:]

        # Normalize observations so the frequency at each round adds up to 1.
        for j in range(self.num_rounds):
            Q = self.y_obs[:,j]*self.degeneracy
            self.y_obs[:,j] = Q/np.sum(Q)

        # Initialize temporary arrays used in objective function
        self.y_calc = np.zeros((self.num_patterns,self.num_rounds),
                               dtype=float)
        self.p = np.zeros((self.num_patterns),dtype=float)

        # Populate initial guesses
        # 1) Try to fit a simple single-site model to the data.  This is the
        #    first guess.
        # 2) If that doesn't converge, assign the initial conc to the frequency
        #    at obs0 and  K to 1.0
        print("Generating initial parameter guesses...",end="")
        sys.stdout.flush()
        self.param_guess = np.zeros((self.num_patterns*2),dtype=float) 
        x_for_guess = np.array(range(rounds_start_at,
                                     rounds_start_at + self.num_rounds),dtype=int)
        for i in range(self.num_patterns):
            y = self.y_obs[i,:]
            conc_guess = self.y_obs[i,0]
            K_guess = 1.0 

            try:
                indep_param, cov = curve_fit(self.indepFit,x_for_guess,y,
                                             p0=(conc_guess,K_guess),maxfev=10000)

                if indep_param[0] <= 0 or indep_param[1] <= 0:
                    raise RuntimeError
                self.param_guess[i] = np.log(indep_param[0]) 
                self.param_guess[self.num_patterns + i] = np.log(indep_param[1])
            except RuntimeError:
                self.param_guess[i] = np.log(1/self.num_patterns) 
                self.param_guess[self.num_patterns + i] = np.log(1.0)

        print("Done.")
        sys.stdout.flush()

        # Create bound list.  conc must be between 0 and 1, K must be positive
        self.bounds = [(None,None) for i in range(self.num_patterns)]   
        self.bounds.extend([(None,None) for i in range(self.num_patterns)])

        #self.bounds[0] = (-1.0,-1.0)
        #self.param_guess[0] = -1.0
        self.bounds[self.num_patterns] = (-1.0,-1.0)
        self.param_guess[self.num_patterns] = -1.0
 
        # Constrain the sum of the initial conc parameters to be 1.
        #self.constraints = ({'type': 'eq',
        #                     'fun': lambda x:  1 - sum(x[:int(len(x)/2)])})
        

    def objective(self,param):
        """
        Objective function to minimize.  Goes like:

        For each pattern i, the observed frequency (theta) at round j is:
 
            theta[i,j] = (theta[i,0]*(beta[i]**j))/sum_over_i_at_j

        Objective function is the rmsd between calculated theta and observed theta
        over all rounds.

        param = [all_thetas...., all_Ks...]
        """

        # Initial state 
        self.p = param[:self.num_patterns]

        # For every round, apply the beta term
        for j, n in enumerate(self.round_exponents):
            if n != 0:
                self.p = self.p*(param[self.num_patterns:])
            self.y_calc[:,j] = self.p/sum(self.p)

        return np.sum(np.power(self.y_calc-self.y_obs,2))


    def objective2(self,param):
        """
        Objective function to minimize.  Goes like:

        ln(theta_x(i)) = A + B*i - ln(sum(Q))
        exp(A + B*i)/sum(Q)

        Objective function is the rmsd between calculated theta and observed theta
        over all rounds.

        param = [all_thetas...., all_Ks...]
        """

        for j, n in enumerate(self.round_exponents):
            self.p = np.exp(self.log_degeneracy + param[:self.num_patterns] + param[self.num_patterns:]*n)
            self.y_calc[:,j] = self.p/np.sum(self.p)

        return np.sum(np.power(self.y_calc-self.y_obs,2))


    def runRegression(self,maxiter=100000000000):
        """
        Run a regression of the objective function given our bound, constraints,
        etc.
        """

        self.start_time = time.time()

        self.fit_result = minimize(fun=self.objective2,x0=self.param_guess,
                                   bounds=self.bounds,
                                   #constraints=self.constraints,
                                   options={"maxiter":maxiter,"maxfun":maxiter})
        self.end_time = time.time()

        print(time.asctime())
        print(self.fit_result.status,self.fit_result.fun,self.fit_result.message)
        print(time.asctime())

    def writeOut(self,out_file):
        """
        Dump the fit result out.
        """

        pickle.dump(self.fit_result,open(out_file,"wb"))

    def returnParam(self):
        """
        """

        initial_theta = self.fit_result.x[:self.num_patterns]
        K_values = self.fit_result.x[self.num_patterns:]

        out = np.zeros((len(K_values),3),dtype=float)
        for i in range(3):
            out[:,i] = np.exp(self.log_degeneracy + initial_theta + K_values*(i+1))
            out[:,i] = out[:,i]/sum(out[:,i])

        return initial_theta, K_values, out   


    def indepFit(self,x,conc=1.0,K=1.0):
    
        return conc*(K**x)
