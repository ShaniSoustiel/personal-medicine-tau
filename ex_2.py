# import libraries
import pandas as pd
import seaborn as sns
import scipy.stats
from statsmodels.formula.api import ols
import statsmodels.api as sm

# input and output paths
input_path = ('/Users/shanishalev/Library/CloudStorage'
              '/OneDrive-Personal/Documents - one drive/Biology TAU'
              '/Year no 03/semester B/personal_medicine_course/hw_2/data/')
output_path = ('/Users/shanishalev/Library/CloudStorage/OneDrive-Personal/'
               'Documents - one drive/Biology TAU/Year no 03/semester B/'
               'personal_medicine_course/hw_2/outputs/')

# load data
gen = pd.read_excel(input_path+'genotypes.xlsx')
phen = pd.read_excel(input_path+'phenotypes.xlsx')

#%% Q1: simple linear regression model in which heterozygous markers are ignored
# formalize the model assumptions
# and calculate P-value scores for th chosen SNP

# filter out the heterozygous markers, cells with "H" in them
gen_no_hetero = gen[gen['check hetrozygoty'] != "Contains H"]





#%% Q2: simple linear regression model in which heterozygous markers are considered
# formalize the model assumptions
# and calculate P-value scores for th chosen SNP

#%% Q3: ANOVA model in which heterozygous markers are ignored
# formalize the model assumptions
# and calculate P-value scores for th chosen SNP
