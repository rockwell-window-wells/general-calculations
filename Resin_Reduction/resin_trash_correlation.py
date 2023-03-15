# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 08:52:30 2023

@author: Ryan.Larson
"""

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.stats.api as sms

df = pd.read_excel("Resin Trash Correlation.xlsx")

# Make cross-tabulation table
ctab = pd.crosstab(df["Reduced or Old Resin Amount"], df["Trash or Repair"])

# Use chi-square test to determine if there is a statistically significant
# relationship between categorical variables
stat, p, dof, expected = stats.chi2_contingency(ctab)

# Interpret p-value
alpha = 0.05
print("Relationship between reduced resin case and trash or repair")
print("p value is " + str(p))
if p <= alpha:
    print("Dependent (reject null hypothesis)")
else:
    print("Independent (do not reject the null hypothesis)")


# Point biserial correlation for going between the measured weight and trash or repair (72s only)
df72 = df[df["Part Size"] == 72]
binary_cols = pd.get_dummies(df72["Trash or Repair"])
df72 = pd.concat((df72, binary_cols["Repair"]), axis=1)

r, p = stats.pointbiserialr(df72["Total Resin Weight"], df72["Repair"])
print("\nCorrelation between measured weight and repairs for 72:")
print("p value is " + str(p))
if p <= alpha:
    print("Dependent (reject null hypothesis)")
else:
    print("Independent (do not reject the null hypothesis)")

# Point biserial correlation between bag number and trash or repair frequency (all sizes)
binary_cols = pd.get_dummies(df["Trash or Repair"])
df = pd.concat((df, binary_cols["Repair"]), axis=1)

r, p = stats.pointbiserialr(df["Bag"], df["Repair"])
print("\nCorrelation between bag numbers and repairs:")
print("p value is " + str(p))
if p <= alpha:
    print("Dependent (reject null hypothesis)")
else:
    print("Independent (do not reject the null hypothesis)")
    
# Mold against repair or trash
ctab = pd.crosstab(df["Mold"], df["Trash or Repair"])
stat, p, dof, expected = stats.chi2_contingency(ctab)
alpha = 0.05
print("\nRelationship between mold and trash or repair")
print("p value is " + str(p))
if p <= alpha:
    print("Dependent (reject null hypothesis)")
else:
    print("Independent (do not reject the null hypothesis)")

