#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 09:57:48 2021

@author: linigodelacruz
"""
# Import packages


# importing all the necessary packages
import numpy as np
import matplotlib.pyplot as plt
  
# importing the style package
from matplotlib import style
  
# creating an array of data for plot
data = np.random.randn(100)
  
# using the style for the plot
#plt.style.use('scientific')
plt.style.use('scientific')

## styles available: ['Solarize_Light2', '_classic_test_patch',
# 'bmh', 'classic', 'dark_background', 'fast', 'fivethirtyeight', 'ggplot', 
#'grayscale', 'scientific', 'seaborn', 'seaborn-bright', 'seaborn-colorblind',
# 'seaborn-dark', 'seaborn-dark-palette', 'seaborn-darkgrid', 'seaborn-deep',
# 'seaborn-muted', 'seaborn-notebook', 'seaborn-paper', 'seaborn-pastel',
# 'seaborn-poster', 'seaborn-talk', 'seaborn-ticks', 'seaborn-white', 
#'seaborn-whitegrid', 'tableau-colorblind10']
  
# creating a plot
plt.plot(data)
  
# show plot
plt.show()

