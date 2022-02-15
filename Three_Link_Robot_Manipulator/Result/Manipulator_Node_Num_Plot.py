import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns; sns.set()
import xlrd
import os
sns.set_style("white")
os.listdir(".")
worksheet = xlrd.open_workbook('SODLRBLS_Manipulator_Node_Num.xls.xls')
sheet_names = worksheet.sheet_names()
c1 = []
for sheet_name in sheet_names:
    sheet = worksheet.sheet_by_name(sheet_name)
    rows = sheet.nrows 
    cols = sheet.ncols 
    all_content = []
    for i in range(0,10):
        cols = sheet.col_values(i)
        c1.append(cols)
fig = plt.figure()
rewards1 = np.array(c1[0][0:25002])
rewards2 = np.array(c1[1][0:25002])
rewards3 = np.array(c1[2][0:25002])
rewards4 = np.array(c1[3][0:25002])
rewards5 = np.array(c1[4][0:25002])
rewards6 = np.array(c1[5][0:25002])
rewards7 = np.array(c1[6][0:25002])
rewards8 = np.array(c1[7][0:25002])
rewards9 = np.array(c1[8][0:25002])
rewards10 = np.array(c1[9][0:25002])
rewards = np.concatenate((rewards1,rewards2,rewards3,rewards4,rewards5,rewards6,rewards7,rewards8,rewards9,rewards10))
episode1 = range(len(rewards1))
episode2 = range(len(rewards2))
episode3 = range(len(rewards3))
episode4 = range(len(rewards4))
episode5 = range(len(rewards5))
episode6 = range(len(rewards6))
episode7 = range(len(rewards7))
episode8 = range(len(rewards8))
episode9 = range(len(rewards9))
episode10 = range(len(rewards10))
episode = np.concatenate((episode1,episode2,episode3,episode4,episode5,episode6,episode7,episode8,episode9,episode10))
sns.lineplot(x=episode,y=rewards)
plt.xlabel("Time(ms)")
plt.ylabel("Number of the ramaining feature neuron")
plt.show()