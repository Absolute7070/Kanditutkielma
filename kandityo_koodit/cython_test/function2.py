#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 08:38:52 2023

@author: degnaiyu
"""
import time 


def function(x, y, n): 
    total = 0 
    for i in range(n): 
        for k in range(n): 
            total += x**y
    return total


start = time.time()
function(19, 12, 100)
end = time.time()



print('time used in python implementation', end - start )
        