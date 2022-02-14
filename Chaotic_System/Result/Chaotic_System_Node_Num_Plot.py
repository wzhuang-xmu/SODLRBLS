import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns; sns.set()
import xlrd
import os
sns.set_style("white")
os.listdir(".")
worksheet = xlrd.open_workbook('SOBLS_Chaotic_System_Node_Num.xls')
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
worksheet = xlrd.open_workbook('SODLRBLS_Chaotic_System_Node_Num.xls')
sheet_names = worksheet.sheet_names()
c2 = []
for sheet_name in sheet_names:
    sheet = worksheet.sheet_by_name(sheet_name)
    rows = sheet.nrows 
    cols = sheet.ncols 
    all_content = []
    for i in range(0,10):
        cols = sheet.col_values(i)
        c2.append(cols)
fig = plt.figure()
rewards1 = np.array(c1[0][0:7002])
rewards2 = np.array(c1[1][0:7002])
rewards3 = np.array(c1[2][0:7002])
rewards4 = np.array(c1[3][0:7002])
rewards5 = np.array(c1[4][0:7002])
rewards6 = np.array(c1[5][0:7002])
rewards7 = np.array(c1[6][0:7002])
rewards8 = np.array(c1[7][0:7002])
rewards9 = np.array(c1[8][0:7002])
rewards10 = np.array(c1[9][0:7002])
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
reward1 = np.array(c2[0][0:7002])
reward2 = np.array(c2[1][0:7002])
reward3 = np.array(c2[2][0:7002])
reward4 = np.array(c2[3][0:7002])
reward5 = np.array(c2[4][0:7002])
reward6 = np.array(c2[5][0:7002])
reward7 = np.array(c2[6][0:7002])
reward8 = np.array(c2[7][0:7002])
reward9 = np.array(c2[8][0:7002])
reward10 = np.array(c2[9][0:7002])
reward = np.concatenate((reward1,reward2,reward3,reward4,reward5,reward6,reward7,reward8,reward9,reward10))
episo1 = range(len(reward1))
episo2 = range(len(reward2))
episo3 = range(len(reward3))
episo4 = range(len(reward4))
episo5 = range(len(reward5))
episo6 = range(len(reward6))
episo7 = range(len(reward7))
episo8 = range(len(reward8))
episo9 = range(len(reward9))
episo10 = range(len(reward10))
episo = np.concatenate((episo1,episo2,episo3,episo4,episo5,episo6,episo7,episo8,episo9,episo10))
sns.lineplot(x=episo,y=reward)
plt.legend(('SOBLS','SODLRBLS'), loc='upper right')
plt.xlabel("Time(ms)")
plt.ylabel("Number of the ramaining feature neuron")
plt.show()