#%% importing things
import sys, os

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from scipy.stats import pearsonr
from scipy.stats import spearmanr

#%% function to read in data
data_dir = r"C:\Users\oyina\Documents\src\measurement_principles\abm-project-\2_Wound infection model editted experiment-spreadsheet.csv"
data = pd.read_csv(data_dir, usecols=range(1,13), skiprows=21)

col_names = data.columns
start_bacteria = [col_name.startswith("count bacteria") for col_name in col_names]
start_neutrophils = [col_name.startswith("count neutrophil") for col_name in col_names]
start_macrophages = [col_name.startswith("count macrophage") for col_name in col_names]
start_patches = [col_name.startswith("count patches") for col_name in col_names]

bacteria_data = data[col_names[start_bacteria]].values.mean(axis=1)
neutrohils_data = data[col_names[start_neutrophils]].values.mean(axis=1)
macrophages_data = data[col_names[start_macrophages]].values.mean(axis=1)
patches_data = data[col_names[start_patches]].values.mean(axis=1)

sum_neutrophils_macrophages = np.sum(neutrohils_data) + np.sum(macrophages_data)
sum_bacteria = np.sum(bacteria_data)
bacteria_measure = sum_bacteria/sum_neutrophils_macrophages

def get_bacteria_measure(data_dir, num_cols):
    data = pd.read_csv(data_dir, usecols=range(1,13), skiprows=21)

    col_names = data.columns
    start_bacteria = [col_name.startswith("count bacteria") for col_name in col_names]
    start_neutrophils = [col_name.startswith("count neutrophil") for col_name in col_names]
    start_macrophages = [col_name.startswith("count macrophage") for col_name in col_names]
    start_patches = [col_name.startswith("count patches") for col_name in col_names]

    bacteria_mean = data[col_names[start_bacteria]].values.mean(axis=1)
    bacteria_std = data[col_names[start_bacteria]].values.std(axis=1)

    neutrohils_mean = data[col_names[start_neutrophils]].values.mean(axis=1)
    neutrohils_std = data[col_names[start_neutrophils]].values.std(axis=1)

    macrophages_mean = data[col_names[start_macrophages]].values.mean(axis=1)
    macrophages_std = data[col_names[start_macrophages]].values.std(axis=1)

    patches_mean = data[col_names[start_patches]].values.mean(axis=1)
    patches_std = data[col_names[start_patches]].values.std(axis=1)

    sum_neutrophils_macrophages_mean = np.sum(neutrohils_mean) + np.sum(macrophages_mean)
    sum_neutrophils_macrophages_std = np.sum(neutrohils_std) + np.sum(macrophages_std)

    sum_bacteria_mean = np.sum(bacteria_mean)
    sum_bacteria_std = np.sum(bacteria_std)

    bacteria_measure = sum_bacteria_mean/sum_neutrophils_macrophages_mean
    # bacteria_measure_std = sum_bacteria_std/sum_neutrophils_macrophages_std
    return bacteria_measure

bact_measure = get_bacteria_measure(data_dir, 13)


#%%  data directories
prob_1_dir = r"C:\Users\oyina\Documents\src\measurement_principles\abm-project-\simulations\prob_1.csv"
prob_8_dir = r"C:\Users\oyina\Documents\src\measurement_principles\abm-project-\simulations\prob_0.8.csv"
prob_6_dir = r"C:\Users\oyina\Documents\src\measurement_principles\abm-project-\simulations\prob_0.6.csv"
prob_4_dir = r"C:\Users\oyina\Documents\src\measurement_principles\abm-project-\simulations\prob_0.4.csv"
prob_2_dir = r"C:\Users\oyina\Documents\src\measurement_principles\abm-project-\simulations\prob_0.2.csv"
treat_prob_1_dir = r"C:\Users\oyina\Documents\src\measurement_principles\abm-project-\simulations\treat_prob_1.csv"
treat_prob_8_dir = r"C:\Users\oyina\Documents\src\measurement_principles\abm-project-\simulations\treat_prob_0.8.csv"
treat_prob_6_dir = r"C:\Users\oyina\Documents\src\measurement_principles\abm-project-\simulations\treat_prob_0.6.csv"
treat_prob_4_dir = r"C:\Users\oyina\Documents\src\measurement_principles\abm-project-\simulations\treat_prob_0.4.csv"
treat_prob_2_dir = r"C:\Users\oyina\Documents\src\measurement_principles\abm-project-\simulations\treat_prob_0.2.csv"

#%% run data through function
bacteria_measure_treat = np.zeros((5))
bacteria_measure = np.zeros((5))

prob_1_measure = get_bacteria_measure(prob_1_dir, 40)
prob_8_measure = get_bacteria_measure(prob_8_dir, 40)
prob_6_measure = get_bacteria_measure(prob_6_dir, 40)
prob_4_measure = get_bacteria_measure(prob_4_dir, 40)
prob_2_measure = get_bacteria_measure(prob_2_dir, 40)
treat_prob_1_measure = get_bacteria_measure(treat_prob_1_dir, 40)
treat_prob_8_measure = get_bacteria_measure(treat_prob_8_dir, 40)
treat_prob_6_measure = get_bacteria_measure(treat_prob_6_dir, 40)
treat_prob_4_measure = get_bacteria_measure(treat_prob_4_dir, 40)
treat_prob_2_measure = get_bacteria_measure(treat_prob_2_dir, 40)

#%% plotting
x_labels = [0.2, 0.4, 0.6, 0.8, 1]
no_treat_measures = [prob_2_measure, prob_4_measure, prob_6_measure, prob_8_measure, prob_1_measure]
treat_measures = [treat_prob_2_measure, treat_prob_4_measure, treat_prob_6_measure, treat_prob_8_measure, treat_prob_1_measure]

fig = plt.figure()
plt.plot(x_labels, no_treat_measures, "o", linestyle="solid", label="without antibiotics")
plt.plot(x_labels, treat_measures, "o", linestyle="solid", label="with antibiotics")
plt.legend()
plt.xlabel("probability of macrophage killing bacteria")
plt.ylabel("bacterial persistence")
# fig.suptitle('test title', fontsize=12)
plt.title("bacterial persistence with changing macrophage killing")
plt.savefig("paper_fig.png")

pearson_treat = pearsonr(x_labels, treat_measures)
pearson_no_treat = pearsonr(x_labels, no_treat_measures)

spearman_treat = spearmanr(x_labels, treat_measures)
spearman_no_treat = spearmanr(x_labels, no_treat_measures)
