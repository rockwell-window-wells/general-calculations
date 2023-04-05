# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 08:52:30 2023

@author: Ryan.Larson
"""

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.stats.api as sms

# Read in data from Excel
df = pd.read_excel("Resin Trash Correlation.xlsx",sheet_name="Sheet1")

#########################
##### DATA CHECKING #####
#########################

# Verify that the resin weights entered are appropriate for the treatment
# specified
resin_old = {36: 33.2, 48: 38.7, 60: 44.2, 72: 51.9, 84: 60.2, 96: 70.7, 102: 80.9}
resin_new = {36: 31.3, 48: 36.5, 60: 41.6, 72: 48.9, 84: 56.7, 96: 66.6, 102: 76.2}

resin_diffs = {}
for key in resin_old.keys():
    resin_diffs[key] = np.around((resin_old[key] - resin_new[key]), 1)
    
df_old = df[df["New or Old Resin Amount"]=="Old"]
df_new = df[df["New or Old Resin Amount"]=="New"]

def get_old_resin(x):
    return resin_old[x]

def get_new_resin(x):
    return resin_new[x]

def get_max_diff(x):
    return resin_diffs[x]

df_old["Expected Resin Weight"] = df_old["Part Size"].apply(get_old_resin)
df_new["Expected Resin Weight"] = df_new["Part Size"].apply(get_new_resin)

df = pd.concat((df_old, df_new), axis=0)
# df.sort_values("Date")
# df.reset_index(inplace=True,drop=True)

df["Diff"] = df["Expected Resin Weight"] - df["Total Resin Weight"]
df["Max Diff"] = df["Part Size"].apply(get_max_diff)

def switch_flag(oldnew,diff,maxdiff):
    if oldnew == "Old":
        if diff > maxdiff/2:
            return "New"
        else:
            return "Old"
    else:
        if diff < -maxdiff/2:
            return "Old"
        else:
            return "New"
        
df["Adjusted New or Old"] = df.apply(lambda x: switch_flag(x["New or Old Resin Amount"], x["Diff"], x["Max Diff"]), axis=1)

# Filter out cases where the part size doesn't remotely match the resin written down
df = df[np.abs(df["Diff"]) <= df["Max Diff"]]
df_old = df[df["Adjusted New or Old"]=="Old"]
df_new = df[df["Adjusted New or Old"]=="New"]
df_old["New Expected Resin Weight"] = df_old["Part Size"].apply(get_old_resin)
df_new["New Expected Resin Weight"] = df_new["Part Size"].apply(get_new_resin)

df = pd.concat((df_old, df_new), axis=0)
df["New Diff"] = df["New Expected Resin Weight"] - df["Total Resin Weight"]

######################
##### STATISTICS #####
######################


# Make cross-tabulation table
ctab = pd.crosstab(df["New or Old Resin Amount"], df["Trash or Repair"])

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
df72 = pd.concat((df72, binary_cols["Good"]), axis=1)
df72 = pd.concat((df72, binary_cols["Repair"]), axis=1)
# df72 = pd.concat((df72, binary_cols["Trash"]), axis=1)

r, p = stats.pointbiserialr(df72["Total Resin Weight"], df72["Good"])
print("\nCorrelation between measured weight and good parts for 72:")
print("p value is " + str(p))
if p <= alpha:
    print("Dependent (reject null hypothesis)")
else:
    print("Independent (do not reject the null hypothesis)")

# r, p = stats.pointbiserialr(df72["Total Resin Weight"], df72["Repair"])
# print("\nCorrelation between measured weight and repair parts for 72:")
# print("p value is " + str(p))
# if p <= alpha:
#     print("Dependent (reject null hypothesis)")
# else:
#     print("Independent (do not reject the null hypothesis)")

# r, p = stats.pointbiserialr(df72["Total Resin Weight"], df72["Trash"])
# print("\nCorrelation between measured weight and trash parts for 72:")
# print("p value is " + str(p))
# if p <= alpha:
#     print("Dependent (reject null hypothesis)")
# else:
#     print("Independent (do not reject the null hypothesis)")


# Point biserial correlation between bag number and trash or repair frequency (all sizes)
binary_cols = pd.get_dummies(df["Trash or Repair"])
df = pd.concat((df, binary_cols["Good"]), axis=1)
df = pd.concat((df, binary_cols["Repair"]), axis=1)
df = pd.concat((df, binary_cols["Trash"]), axis=1)

r, p = stats.pointbiserialr(df["Bag"], df["Good"])
print("\nCorrelation between bag numbers and good parts:")
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

print("\n\n")