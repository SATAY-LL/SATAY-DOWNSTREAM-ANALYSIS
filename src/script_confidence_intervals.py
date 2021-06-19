# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 17:09:36 2021

@author: floor
"""


from src.python_modules.module_confidence_intervals import confidenceintervals

from tkinter import Tk
from tkinter import filedialog

        
        



#%% Selecting datasets from your local files 

root = Tk()
root.filename =  filedialog.askopenfilename(title = "pergene insertions file replica 1",filetypes = (("txt files","*.txt"),("all files","*.*")))
filename_replica_a=root.filename
root.withdraw()        

root = Tk()
root.filename =  filedialog.askopenfilename(title = "pergene insertions file replica 2",filetypes = (("txt files","*.txt"),("all files","*.*")))
filename_replica_b=root.filename
root.withdraw() 
        
        
test = confidenceintervals(filename_replica_a, filename_replica_b)
   
        
        