#!/usr/bin/env python
# coding: utf-8

# # Auxiliary functions 

# In[12]:



import re 



# ## Read file to a string 

# In[13]:


def readFileToString(absolutePath: str) -> str:  
    '''
    Read file by absolute path. Returns string 
    
    
    Parameters
    ----------
        absolutePath: absolute path of the file 
        
    
    Return
    ------
        Returns the content of the file as string 
    
    '''
    
    file = open(absolutePath, 'r') 
    filecontent = file.read()   
    file.close()     # closing the file 
    return filecontent


# ## Finding variables

# In[14]:


def findVariables(exprStr: str) -> tuple:
    '''
    Extract variables in expression string. 
    
    Parameters
    ---------
        exprStr: expression as string. Supposed to be Python-compatible formula (at least Fortran, but not tested)
        
    
    Return
    ------
    Tuple of two lists. 
    
        First one: alphanumerically-ordered coordinate list e.g.
        ['x1', 'x2', 'x3', 'y1', 'y2', 'y3', 'z1', 'z2', 'z3']
        
        Second one: alphanumerically-ordered  parameters list 
        e.g. ['A1', 'A2']
        
    '''
    
    patternForCoordinates = r'[xyz]{1}[123]{1}'                   # pattern for: x1, x2, y1 etc. 
    patternForParameters = r'[A-Z]{1}\d{1}'                # pattern for: A1, A2, ... Z9  
    
    coordinateList =sorted(list(set( re.findall(patternForCoordinates, exprStr)  )))
    parametersList =sorted(list(set(  re.findall(patternForParameters, exprStr) )))
    
    return coordinateList, parametersList 
    
    


# ## Regrouping coordinates 

# Regrouping the coordinates from the output of the previous function `findVariables` into form [(x1, y1, z1), (x2, y2, z2)]: 

# In[2]:


def regroupToCoordinateTriple_findVariables(coordinateList: list, isCoordinateTripletInTuple: bool = False ) -> list: 
    '''
    Regrouping the coordinates from the output of the previous function `findVariables` into form [(x1, y1, z1), (x2, y2, z2)]
    
    Parameters
    ----------
        coordinateList: output coordinates from findVariables-function 
        
        isCoordinateTripletInTuple: whether the output coordinate triples be separated by parentheses (see Return)
        
    Return
    ------
        list of traditional coordinatetriplets [(x1, y1, z1), (x2, y2, z2)] or in ordered list without parentheses [x1, y1, z1, x2, y2, z2 ]
        
    Warning
    -------
        User check that the coordinateList is the output of the function `findVariables` 
    '''
    
    tripletList = []   # save all triplets 
    
    
    # each particle's components 
    xCoordinates = []
    yCoordinates = []
    zCoordinates = []
    
    
    # separating each particle's component into their lists 
    for component in coordinateList: 
        if 'x' in component: 
            xCoordinates.append(component) 
        elif 'y' in component: 
            yCoordinates.append(component)
        elif 'z' in component:   
            zCoordinates.append(component)
        else: 
            raise Exception("No other than x,y,z components")
    
    
    
    
    if isCoordinateTripletInTuple: 
        # forming triplets 
        for x,y,z in zip(xCoordinates, yCoordinates, zCoordinates): 
            tripletList.append((x,y,z))
    else: 
        # forming ordered list 
        for x,y,z in zip(xCoordinates, yCoordinates, zCoordinates): 
            tripletList.append(x)
            tripletList.append(y)
            tripletList.append(z)
        
        
    return tripletList
    
    
    
    
    


# ## Error handling 

# ### check number of initial guess is same as number of parameters found 

# In[3]:


def check_numberOfIni_VS_numberOfFound_parameters(initialGuessList:list, foundParameters: list): 
    if len(initialGuessList) != len(foundParameters): 
        raise Exception("Number of initial guess not matching number of found parameters")


# In[ ]:




