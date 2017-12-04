import math
import random
import logging
import os, shutil, sys
import numpy as np

# get initial coordinates of xsf file
def get_data(filename):
    with open(filename,'r') as fin:
        lines=fin.readlines()
    tot_num=int(lines[8].split()[0])
    energy=float(lines[0].split()[4])
    cell=np.zeros((3,3))
    #get box dimension
    for i in range(3):
        for j in range(3):
            cell[i][j]=lines[4+i].split()[j]
    #read atom coordinates, when 'Al' +1 to count when 'Ti' appear recount Ti type atom num,store x y z and fx fy fz if possible
    cor=np.zeros((tot_num,6))
    count=0
    atom_type=1
    atom_type_num=[]
    a=str(lines[9].split()[0])
    for i in range(tot_num):
        for j in range(1,7):
            cor[i][j-1]=lines[9+i].split()[j]
        if (a==str(lines[9+i].split()[0])):count+=1
        else :
        #elif (a<>str(lines[9+i].split[0])):
            atom_type_num.append(count)
            a=lines[9+i].split()[0]
            count=1
            atom_type+=1
    atom_type_num.append(count)
    return cell,cor

def eval_func(filename):
    cell=get_data(filename)[0]
    cor=get_data(filename)[1]

    score = 1
    f = open('conc.in', 'w')
    f.write('1.0 1.0 1.0 90 90 90' + '\n')
    f.write('-5.0 5.0 0' + '\n')
    f.write('-5.0 -5.0 10.0' + '\n')
    f.write('8.0 8.0 8.0' + '\n')

    for i in range(4320):
        spec = 'Al'
        f.write(str(cor[i][0]/cell[0][0]*10.0) + '  ' +str(cor[i][1]/cell[1][1]*10.0) + '  '+str(cor[i][2]/cell[2][2]*8.0) + '  '+ spec + '\n')

    for i in range(4320,4800):
        spec = 'Ti'
        f.write(str(cor[i][0]/cell[0][0]*10.0) + '  ' +str(cor[i][1]/cell[1][1]*10.0) + '  '+str(cor[i][2]/cell[2][2]*8.0) + '  '+ spec + '\n')
    f.close()

    comline1 = '/home/mingfei/bin/corrdump -l lat.in -2=3.0 -3=1.25 -clus'
    comline3 = '/home/mingfei/bin/corrdump -noe -2=3.0 -3=1.25 -s=conc.in > tcorr_final.out'



    os.system(comline1)
    os.system(comline3)

    f = open('target_correlations', 'r')
    randcorr = f.readlines()
    f.close()
    
    v1 = np.zeros((len(randcorr[0].strip().split()),1))
    c1 = 0
    for i in randcorr[0].strip().split():
        v1[c1, 0] = np.float(i)
        c1 = c1 + 1



    f = open('tcorr_final.out', 'r')
    trialcorr = f.readlines()
    f.close()

#   v1 = np.zeros((len(trialcorr[0].strip().split()),1))
    v2 = np.zeros((len(trialcorr[0].strip().split()),1))

#   c1 = 0
#   for i in randcorr[0].strip().split():
#       v1[c1, 0] = np.float(i)
#       c1 = c1 + 1

    c1 = 0
    for i in trialcorr[0].strip().split():
        v2[c1, 0] = np.float(i)
        c1 = c1 + 1      

#    print v2

    CR1 = np.zeros((200,1))
    CD1 = np.zeros((200,1))
    CN1 = np.zeros((200,1))
    CO1 = np.zeros((200,1))

    f = open('clusters.out','r')
    c1 = f.readlines()
    nl = sum(1 for line in open('clusters.out'))
    f.close()

    i=0
    j=0
    while (i< nl):
        CD1[j] = int(c1[i])
        R = float(c1[i+1])
        CR1[j] = R
        if R == 0.0:
            CO1[j] = 0.0
        else:
            CO1[j] = 1/R**6
        N = int(c1[i+2])
        i = i + 4 + N
        CN1[j] = N
        j = j+1

    CR = CR1[0:j]
    CD = CD1[0:j]
    CN = CN1[0:j]
    CO = CO1[0:j]

    DV = np.zeros((j,1))

    for i in range(0,j):
        DV[i] = (float(v2[i])-float(v1[i]))*float(CD[i])*float(CO[i])
       
    #score = np.linalg.norm(v2-v1)
    score = np.linalg.norm(DV)    
    return float(score)
for i in range(25,51):
    filename='out_'+str(i*500).zfill(6)+'.xsf'
    a=eval_func(filename)
    print str(i*500).zfill(6),a

