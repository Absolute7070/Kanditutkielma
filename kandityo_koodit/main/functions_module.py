#!/usr/bin/env python
# coding: utf-8

# # Metropolis, Integration, gradient descent 

# In[ ]:


import numpy as np 
import matplotlib.pyplot as plt 


import sympy as sym 
from sympy import symbols, diff, lambdify 


import parameters as pr  

from numba import jit



import warnings 
warnings.filterwarnings("error")



# ![tensori_yksittainen.png](attachment:d431eaab-739f-4fe2-b30b-c81d6997d46f.png)

# ## Metropolis-sampling 

# ### Thermalisation 

# In[3]:


# testattu 

def metropolisSamplingFunction(
    coordinateValueRange: list or tuple, numberOfParticles: int, 
    numberOfConfig: int, 
    probabilityFunction, params: list, 
    numberOfIterations: int = 100, 
    delta: int  = 0.5) -> np.ndarray: 
    
    '''
    Generating samples according to Metropolis-algorithm for at least 2 particles. Does not work for 1 particle. 
    
    Parameters
    ----------
        coordinateValueRange: lower limit and upper limit for generated coordinates. E.g. (-10, 10) in units meters 
        
        numberOfParticles: how many particles 
        
        numberOfConfig: how many possible configurations/samples for each particle 
        
        numberOfIterations: how many iterations of the whole configuration space 
        
        delta: determining the step size for each particle's coordinates
    

    
    Return 
    ------
    
        tensor with shape (size, number of samples/configurations, 3 coordinates) (See picture above)
    '''
    
    # see drawing in a sketch in kandi folder of the tensor 
    initial_configSamples =  np.random.uniform(coordinateValueRange[0], coordinateValueRange[1], size = (numberOfParticles, numberOfConfig, 3)) # e.g. size is 2 particle, 10 configurations/samples, 3 coordinates for each particle 
    
    
    # how many iterations over the whole configuration space 
    for iteration in range(numberOfIterations): 
    
        for configIndex in range(initial_configSamples.shape[1]): 
             #  current configurations of each particle in format: (each row representing different particle's coordinates)
           #  array([[-0.5583685 , -0.04608995,  0.15500853],
           # [ 0.66255653,  0.66301583, -0.85159876]])
            currentConfig = initial_configSamples[:, configIndex, :]               

            # move particle individually: when suggesting a move for one particle, the other particles remain unmoved. 
            for index, currentParticleCoordinates in zip( list(range(len(currentConfig))), currentConfig): 

                deltai = np.random.uniform(-delta, delta, 3)    # generate step size for each coordinate separately 

                coordinates_trial = currentParticleCoordinates + deltai 

                # defining trial configuration where the coordinate triplet for the current particle is replaced by trial coordinates 
                trialConfiguration = currentConfig.copy() # make copy of original configuration 
                trialConfiguration[index, : ] = coordinates_trial # this is trial configuration with the current position replaced by trial position
                
                
               
                try: 
                    w = probabilityFunction(  trialConfiguration, params  )/probabilityFunction(currentConfig, params)
                except RuntimeWarning: 
                    print('Virhe tapahtui w laskemisessa')
                    return 
                
                
                if w >= 1: 
                    currentConfig[index] = coordinates_trial  # accept the trial for the current particle 
                elif w < 1: 
                    r = np.random.random()   # generating a random number between [0, 1)  
                    if r <= w: 
                        currentConfig[index] = coordinates_trial
                # else the current particle's coordinates are not replaced 


            # update this configuration set 
            initial_configSamples[:, configIndex, :] = currentConfig
    
    return initial_configSamples


# ### Step suggestion 

# In[ ]:


def metropolisStepSuggestion(configuration: np.ndarray, 
                             nthParticle: int, 
                             probabilityFunction, 
                             params: list,  
                            delta: int  = 0.5) -> np.ndarray: 
    
    '''
    Making single suggestion on one particle in n-th position
    
    Parameters
    ----------
        configuration: current configuration of particles in form
                        array([[-0.5583685 , -0.04608995,  0.15500853],
                           [ 0.66255653,  0.66301583, -0.85159876]]). 
                        Each row representing one particle's position 
                        
        nthPosition:
                    To which particle in the table (configuration) we are making 
                    metropolis displacement suggestion.
                    E.g. making suggestion on 1. particle, nthPosition = 0, since its position
                    is on the first row of the configuration table. 
        
        probabilityFunction: 
                    probability function defined above 
        
        params: 
                list of parameters as float 
                
        delta: 
                step size
                
    Return
    -------
        The original configuration table with n-th particle having gone through the Metropolis-suggestion
        (displaced or not)
        
            
    '''
    
    # get the particle of interest 
    particleOfInterest = configuration[nthParticle]
    
    # generating step size 
    deltai = np.random.uniform(-delta, delta, 3) 
    
    # trial position 
    coordinates_trial = particleOfInterest + deltai 
    
    # trial configuration 
    trialConfiguration = configuration.copy() # make copy of original configuration 
    trialConfiguration[nthParticle] = coordinates_trial  
    
    
    # compute acceptability 
    w = probabilityFunction(  trialConfiguration, params  )/probabilityFunction(configuration, params)
    
    # Metropolis-test 
    if w >=1: 
        return trialConfiguration 
    elif w <1: 
        r = np.random.random()   # generating a random number between [0, 1)  
        if r <= w: 
            return trialConfiguration 
    
    # if change not accpeted, then return original configuration table 
    return configuration 

    
    
    


# ## Monte Carlo energy integral approximation 

# In[3]:


# testattu 
def monteCarloIntegrationFunction( samples, localEnergyFunction, params: list, needError: bool = False): 
    
    '''
    
    Compute ground state energy approximation. Based on samples from previous functions. 
    
    
    Parameters
    ---------
    
        samples: tensor from the return value of metropolisSamplingFunction
     
        localEnergyFunction: energy function 
        
        params: list of parameters, either symbolic or as values 
        
        needError: need for computing error. If symbolic parameters, set False. 
    
    
    
    Return
    -------
    tuple or float 
    
        if needError=True: returns 3-tuple with (energy approximation, variance std error of mean)
        if needError=False: return energy approximation only as integer 
    
        if symbolic parameters: returns sympy expression with symbolic parameters 
        if parameters are as values: return energy 
    
    '''
    
    
    N = samples.shape[1]    # how many different configurations 
    
     # in the sum: for each configuration of the electrons (i.e. 2d array, each row representing one electron's positions), we apply the local energy function
    energyApprox = 1/N * np.sum(     np.array( [localEnergyFunction(samples[:, i, :], params) for i in range(samples.shape[1])] )     ) 
    
    
    if needError == True: 
        #variance 
        variance =  1/N * np.sum(     np.array( [localEnergyFunction(samples[:, i, :], params) for i in range(samples.shape[1])] ) **2     ) - energyApprox**2 # < f^2 >- < f >^2 

        # standard error of means 
        stdErrorOfMeans = np.sqrt(variance/N)
        
        
        return energyApprox, variance, stdErrorOfMeans
    
    else: 
        return energyApprox
    
    
   
    
    

    


# ## Gradient Descent 

# In[4]:


# testattu 

def gradientDescentFunc(
    expr, 
    paramVariables: list, 
    initialParamValues: list, 
    learning_rate: float = 0.001,  
    limitForLoop: int = 100, 
    energyGradientScaleToBe: float = 0.001): 
    '''
    Compute gradient descent from sympy expression, which has symbols. 
    
    Gradient descent update algorithm: alpha_i = alpha_(i-1)- learning_rate * gradient(expr)
    
    - expr: sympy expression 
    - paramVariables: list of sympy symbols used in the "expr". e.g. [sym.Symbol('a'), sym.Symbol('b')]
    - paramValues: current values of the parameters. Acts as initial values.
    - learning rate 
    - limitForLoop: maximum iterations 
    - energyGradientScaleToBe: the requirement you put on energy gradient, s.t it is near zero optimally. For zero gradient, the parameter is fully minimized and converged. 
    
    
    return: 
         list of entered optimized parameters as float 
    
    Example: simple usage 
    
        # Define the sympy symbols to be used in the function
        x = symbols('x')
        y = symbols('y')
        #Define the function in terms of x and y
        f1 = (x-2) ** 2 + (y-2)**2+5      # sympy expression 


        paramVariables = [x, y]
        initialParamValues = [3.0, 3.0]

        gradientDescentFunc(f1, paramVariables = paramVariables, initialParamValues = initialParamValues )
    
    
    '''
    
    
    # error check: whether the length of paramVariables is the same as initialParamValues 
    if (len(paramVariables) != len(initialParamValues)): 
        print('Error: Length of the list of initial values must be the same as number of symbolic variables!!!')
        return 
    
    
    partialDerivativeList = [diff(expr, variable) for variable in paramVariables ]   # list of partial derivatives with respect to possible variables in the expr
    
    paramValuesList = initialParamValues.copy()  # list of parameter values, each representing its own parameter e.g. [2.0, 2.0] for [alpha, beta]
    
    
    # dictionary with format: 
    #     {
    #         sym.Symbol('a'):3, 
    #         sym.Symbol('b'):4, 

    #     }
    
    
    
    
    # list of gradients for determining when to cut the loop 
    # if parameters reached minimum, gradient should be zero
    gradientList = np.random.randint(low = 100, size = len(paramValuesList))  # initiating gradientList with random large numbers

    
    
    counter = 0 
    #Perform gradient descent
    while ( np.array(gradientList) < energyGradientScaleToBe).all() == False :         # if all gradients are bigger or equal to 0.1, continue descending 
        counter += 1 
        
        # clean the list
        gradientList = []
        
        # update each parameter value with negative gradient descent 
        for index in range(len(paramValuesList)): 
            # for later substitution when evaluating derivatives  
            substitutionDict = {paramVariables[i]:paramValuesList[i] for i in range(len(paramVariables))}
                
                
            # gradient-descenting one parameter 
            paramValuesList[index] -= partialDerivativeList[index].evalf(subs=substitutionDict)*learning_rate
            
            gradientList.append(partialDerivativeList[index].evalf(subs=substitutionDict))  # appending the gradient to the list 
        

        
        # if user not wanting for too much loops, or not wanting to reach the local minimum,
        # or loop last too long, break the loop 
        if counter >= limitForLoop: 
            break 
    
    
    return paramValuesList

    


# In[ ]:




