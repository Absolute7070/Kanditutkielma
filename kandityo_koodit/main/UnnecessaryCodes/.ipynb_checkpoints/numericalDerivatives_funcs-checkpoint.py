#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np 


# # Gradient

# In[ ]:


def gradientNumFunc(potentialFunc) -> tuple: 
    '''
    Use np.gradient to compute gradient numerically. 
    
    Args: 
        - potentialFunc: the function to which we are taking the gradient. 
    returns: 
        tuple in with three components. e.g (Ex,Ey,Ez)
    
    
    Example: 

        x,y,z = np.mgrid[0:20:1, 0:20:1, 0:20:1]

        V = 2*x**2 + 3*y**2 - 4*z # just a random function for the potential

        Ex,Ey,Ez = np.gradient(V)
        
        Ex[0,0,0], Ey[0,0,0], Ez[0,0,0]     # referring compute the components for x=0, y=0, z=0
    
    '''
    
    return np.gradient(potentialFunc)


# # Laplacian 

# In[4]:


def laplacianNumFunc(potentialFunc) -> np.ndarray: 
    '''
     Use np.gradient to compute laplacian numerically. 
     
     
     Args: 
        - potentialFunc: the function to which we are taking the gradient. 
     returns: 
        tuple in with three components. e.g (Ex,Ey,Ez)
        
        
        
    example: 
        x,y,z = np.mgrid[0:20:1, 0:20:1, 0:20:1]

        V = 2*x**2 + 3*y**2 - 4*z # just a random function for the potential

        Ex,Ey,Ez = np.gradient(V)

        laplacian_V = np.gradient(Ex)[0]+ np.gradient(Ey)[1] + np.gradient(Ez)[2]

        
    '''
    Ex,Ey,Ez = np.gradient(potentialFunc) 
    laplacianOfpotentialFunc = np.gradient(Ex)[0]+ np.gradient(Ey)[1] + np.gradient(Ez)[2]
    
    return laplacianOfpotentialFunc
    


# In[ ]:





# In[ ]:




