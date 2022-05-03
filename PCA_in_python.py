# -*- coding: utf-8 -*-
"""
Created on Mon May  2 14:43:22 2022

@author: Ryan.Larson
"""

# Code modified from tutorial at https://towardsdatascience.com/pca-clearly-explained-how-when-why-to-use-it-and-feature-importance-a-guide-in-python-7c274582c37e#:~:text=PCA%20technique%20is%20particularly%20useful,for%20denoising%20and%20data%20compression.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from sklearn import datasets
from sklearn.decomposition import PCA
import pandas as pd
from sklearn.preprocessing import StandardScaler

def biplot(score, coeff , y):
    '''
    Author: Serafeim Loukas, serafeim.loukas@epfl.ch
    Inputs:
       score: the projected data
       coeff: the eigenvectors (PCs)
       y: the class labels
   '''
    xs = score[:,0] # projection on PC1
    ys = score[:,1] # projection on PC2
    n = coeff.shape[0] # number of variables
    plt.figure(figsize=(10,8), dpi=100)
    classes = np.unique(y)
    colors = ['g','r','y']
    markers=['o','^','x']
    for s,l in enumerate(classes):
        plt.scatter(xs[y==l],ys[y==l], c = colors[s], marker=markers[s]) # color based on group
    for i in range(n):
        #plot as arrows the variable scores (each variable has a score for PC1 and one for PC2)
        plt.arrow(0, 0, coeff[i,0], coeff[i,1], color = 'k', alpha = 0.9,linestyle = '-',linewidth = 1.5, overhang=0.2)
        plt.text(coeff[i,0]* 1.15, coeff[i,1] * 1.15, "Var"+str(i+1), color = 'k', ha = 'center', va = 'center',fontsize=10)

    plt.xlabel("PC{}".format(1), size=14)
    plt.ylabel("PC{}".format(2), size=14)
    limx= int(xs.max()) + 1
    limy= int(ys.max()) + 1
    plt.xlim([-limx,limx])
    plt.ylim([-limy,limy])
    plt.grid()
    plt.tick_params(axis='both', which='both', labelsize=14)
    
    

def main():
    plt.style.use('ggplot')
    
    # Load the data
    iris = datasets.load_iris()
    X = iris.data
    y = iris.target
    
    # Z-score the features
    scaler = StandardScaler()
    scaler.fit(X)
    X = scaler.transform(X)
    
    # The PCA model
    pca = PCA(n_components=2) # estimate only 2 PCs
    X_new = pca.fit_transform(X) # project the original data into the PCA space
    
    # Plot the data before and after the PCA transform and also color code each
    # point (sample) using the corresponding class of the flower (y).
    fig, axes = plt.subplots(1,2)
    
    axes[0].scatter(X[:,0], X[:,1], c=y)
    axes[0].set_xlabel('x1')
    axes[0].set_ylabel('x2')
    axes[0].set_title('Before PCA')
    
    axes[1].scatter(X_new[:,0], X_new[:,1], c=y)
    axes[1].set_xlabel('PC1')
    axes[1].set_ylabel('PC2')
    axes[1].set_title('After PCA')
    
    # plt.show()
    
    print("\nExplained variance (eigenvalues): {}".format(pca.explained_variance_))
    
    print("\nExplained variance (% total variance): {}".format(pca.explained_variance_ratio_))
    
    # Number of rows in pca.components_ is the number of principal components.
    # Number of columns in pca.components_ is the number of features
    print("\nPrincipal components: \n{}".format(abs(pca.components_)))
    # The result is below:
    # [[ 0.52106591  0.26934744  0.5804131   0.56485654]
    #  [ 0.37741762  0.92329566  0.02449161  0.06694199]]
    # Looking at PC1:
    # [ 0.52106591  0.26934744  0.5804131   0.56485654]
    # Features 1, 3, and 4 are most important for PC1
    
    mpl.rcParams.update(mpl.rcParamsDefault) #reset ggplot style
    
    # Call the biplot function for only the first 2 PCs
    biplot(X_new[:, 0:2], np.transpose(pca.components_[0:2, :]), y)
    plt.show()
    
    
if __name__ == "__main__":
    main()