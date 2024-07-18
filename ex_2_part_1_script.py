# import libraries
import pandas as pd
import seaborn as sns
from statsmodels.formula.api import ols
import statsmodels.api as sm
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

#%% define input and output paths

input_path = ('/Users/shanishalev/Library/CloudStorage'
              '/OneDrive-Personal/Documents - one drive/Biology TAU'
              '/Year no 03/semester B/personal_medicine_course/hw_2/data/')
output_path = ('/Users/shanishalev/Library/CloudStorage/OneDrive-Personal/'
               'Documents - one drive/Biology TAU/Year no 03/semester B/'
               'personal_medicine_course/hw_2/outputs/')

#%% load data

gen = pd.read_excel(input_path+'genotypes.xlsx')
phen = pd.read_excel(input_path+'phenotypes_1.xlsx')
phen_name = phen.iloc[787, 1]
SNP =

#%%
# Select row 780 and columns 7-99
selected_row = phen.iloc[787, 7:100]

# Drop NaN values and convert to list
values_list = selected_row.dropna().tolist()

print(values_list)

#%% Q1: simple linear regression model in which heterozygous markers are ignored
# formalize the model assumptions
# and calculate P-value scores for th chosen SNP

# filter out the heterozygous markers, cells with "H" in them
# gen_no_hetero = gen[gen['check hetrozygoty'] != "Contains H"]



# Phenotype data (Bmax values in pmol/mg)
pheno = np.array(values_list)

# Genotype data ('B' and 'D')
genes = np.array(['B', 'D', 'B', 'B', 'D', 'D', 'B', 'D', 'B', 'B',
                  'B', 'D', 'B', 'B', 'B', 'D', 'B', 'D', 'B', 'B',
                  'D', 'B', 'B'])

# Remove missing genotype values (NaN)
X = np.where(genes == 'B', 0, np.where(genes == 'D', 1, np.nan))
Y = pheno[~np.isnan(X)]

# Calculate regression parameters
Xav = np.mean(X)
Yav = np.mean(Y)
Sxx = np.sum(X**2) - len(X) * Xav**2
Sxy = np.dot(X, Y) - len(X) * Xav * Yav
b1 = Sxy / Sxx
b0 = Yav - b1 * Xav
Yhat = X * b1 + b0

# Calculate F-statistic
n = len(X)
k = 2
SST = np.sum((Y - Yav)**2)
SSE = np.sum((Y - Yhat)**2)
SSR = np.sum((Yhat - Yav)**2)
MST = SST / (n - 1)
MSE = SSE / (n - k)
MSR = SSR / (k - 1)
F = MSR / MSE

# Calculate degrees of freedom
df1 = k - 1
df2 = n - k

# Calculate p-value
alpha = 0.05
p_value = 1 - stats.f.cdf(F, df1, df2)

# Print results
print(f"R-squared (ro2): {1 - SSE / SST:.4f}")
print(f"F-statistic: {F:.4f}")
print(f"Degrees of freedom (numerator, denominator): ({df1}, {df2})")
print(f"P-value: {p_value:.6f}")

# Plot regression line
plt.figure(figsize=(8, 6))
plt.scatter(X, Y, label='Data points', color='b')
plt.plot(X, Yhat, label='Regression line', color='r')
plt.xlabel('Genotype value [B=0, D=1]')
plt.ylabel('Phenotype value [Bmax, pmol/mg]')
plt.title('Linear Regression without Heterozygous')
plt.grid(True)
plt.legend()
plt.show()


#%%

model_0 = ols(formula="sepal_width ~ sepal_length", data=iris_data).fit()
intercept, slope = model_0.params
predictions = intercept + slope*iris_data["sepal_length"]
plt.scatter(iris_data["sepal_length"],iris_data["sepal_width"])
plt.plot(iris_data["sepal_length"],predictions, color='red')


#%% Q2: simple linear regression model in which heterozygous markers are considered
# formalize the model assumptions
# and calculate P-value scores for th chosen SNP

#%% Q3: ANOVA model in which heterozygous markers are ignored
# formalize the model assumptions
# and calculate P-value scores for th chosen SNP
