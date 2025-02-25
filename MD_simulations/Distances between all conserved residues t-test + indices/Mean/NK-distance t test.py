# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 17:04:25 2023

@author: felderfl
"""
import matplotlib.pyplot as plt
from scipy import stats

#means of the NK distances of all the MD data

defluorinating_data = [1.9751216940631227, 1.8835529157947675, 1.9764182177610998, 1.9873374488387192, 2.0977577316200855, 2.193737072108241, 2.1282794968034837, 2.2185140596537027, 2.125891187641151, 2.0021879314718167, 2.03989717434993]
nondefluorinating_data = [1.9827825361277414, 1.9949400181216628, 1.9569897025849592, 1.9949400181216628, 2.125891187641151, 2.153160183760938, 2.293730880423418]

# Perform an independent two-sample t-test
t_stat, p_value = stats.ttest_ind(defluorinating_data, nondefluorinating_data)

# Print the results
print("t-statistic:", t_stat)
print("p-value:", p_value)

# Boxplot
data = [defluorinating_data, nondefluorinating_data]
labels = ['Defluorinating (n={})'.format(len(defluorinating_data)), 'Non-defluorinating (n={})'.format(len(nondefluorinating_data))]

plt.figure(figsize=(10, 6))
plt.boxplot(data, labels=labels)
plt.title('Boxplot of Defluorinating vs Non-defluorinating Data\np-value: {:.4f}'.format(p_value))
plt.ylabel('Values')
plt.xlabel('Groups')


# Increase font size of sample sizes
plt.xticks(fontsize=12)

plt.show()