# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 09:47:45 2023

@author: royno
"""

# load KEGG
import scipy.io
import pandas as pd
filename = r'X:\Common\useful_datasets\kegg_pathways.mat'

mat = scipy.io.loadmat(filename,simplify_cells=True)

mat = mat["kegg_pathways"]


df = pd.DataFrame({"pathway":[""]*len(mat),"genes":[""]*len(mat)})
for i in range(len(mat)):
    df["pathway"][i] = mat[i]["name"][5:]    
    df["genes"][i] = "_".join(mat[i]["genes"])
    
df.to_csv(r'X:\Common\useful_datasets\kegg_pathways.csv',index=False)