#!/usr/bin/python
import numpy as np
import random as rd
import os
#import gen
rd.seed()

def calc_dist(lst,n1,n2,H):
    if  (lst[n2][0]-lst[n1][0]>=(0.5*H[0][0])):
        a=(lst[n1][0]-lst[n2][0]+H[0][0])**2
    elif  lst[n2][0]-lst[n1][0]<(-0.5*H[0][0]):
        a=(lst[n1][0]-lst[n2][0]-H[0][0])**2
    else:
        a=(lst[n1][0]-lst[n2][0])**2
    if  lst[n2][1]-lst[n1][1]>=(0.5*H[1][1]):
        b=(lst[n1][1]-lst[n2][1]+H[1][1])**2
    elif  lst[n2][1]-lst[n1][1]<(-0.5*H[1][1]):
        b=(lst[n1][1]-lst[n2][1]-H[1][1])**2
    else:
        b=(lst[n1][1]-lst[n2][1])**2
    if  lst[n2][2]-lst[n1][2]>=0.5*H[2][2]:
        c=(lst[n1][2]-lst[n2][2]+H[2][2])**2
    elif  lst[n2][2]-lst[n1][2]<(-0.5*H[2][2]):
        c=(lst[n1][2]-lst[n2][2]-H[2][2])**2
    else:
        c=(lst[n1][2]-lst[n2][2])**2
    return np.sqrt(a+b+c)


class MonteCarlo(object):
# get initial coordinates of xsf file
    def get_data(self):
        with open(self.filename,'r') as fin:
            lines=fin.readlines()
        self.tot_num=int(lines[8].split()[0])
        self.energy=float(lines[0].split()[4])
        #get box dimension
        for i in range(3):
            for j in range(3):
                self.cell[i][j]=lines[4+i].split()[j]
        #read atom coordinates, when 'Al' +1 to count when 'Ti' appear recount Ti type atom num,store x y z and fx fy fz if possible
        cor=np.zeros((self.tot_num,6))
        count=0
        atom_type=1
        atom_type_num=[]
        a=str(lines[9].split()[0])
        for i in range(self.tot_num):
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
        self.data = cor
        del cor
        self.atom_type=atom_type
        self.atom_type_num=map(int,atom_type_num)
        return
    #print output xsf file for ANN
    def xsf_new(self):
        ii=str(self.tmp_name)
        energy=self.energy
        file_out=str(ii)
        fout=open(file_out,'w')
        saaa=str('# total energy = '+ str(energy) +' eV\n\n')
        fout.write(saaa)
        fout.write('CRYSTAL\nPRIMVEC\n')
        for i in range(3):
            fout.write(str(str(self.cell[i][0])+' '+str(self.cell[i][1])+' '+str(self.cell[i][2])+'\n'))
        fout.write(str('PRIMCOORD\n'+str(self.tot_num)+' 1\n'))
        for i in range(int(self.atom_type_num[0])):
            fout.write(str('Al  '+str("%0.4f" % self.data[i][0])+'   '+str("%0.4f" % self.data[i][1])+'   '+str("%0.4f" % self.data[i][2])+'   '+str(self.data[i][3])+'   '+str(self.data[i][4])+'   '+str(self.data[i][5])+'\n'))
        for i in range(int(self.atom_type_num[0]),self.tot_num):
            fout.write(str('Ti  '+str("%0.4f" % self.data[i][0])+'   '+str("%0.4f" % self.data[i][1])+'   '+str("%0.4f" % self.data[i][2])+'   '+str(self.data[i][3])+'   '+str(self.data[i][4])+'   '+str(self.data[i][5])+'\n'))
        return
 
    def lmp_new(self):
        with open('new.lmp','w') as fout:
            fout.write('# elastic\n\n')
            fout.write('    '+str(self.tot_num)+'    atoms\n     2 atom types\n\n')
            fout.write('    0.0     '+str(self.cell[0][0])+' xlo xhi\n')
            fout.write('    0.0     '+str(self.cell[1][1])+' ylo yhi\n')
            fout.write('    0.0     '+str(self.cell[2][2])+' zlo zhi\n\n')
            fout.write(' %0.4f     %0.4f     %0.4f   xy xz yz\n\n'%(self.cell[0][1],self.cell[0][2],self.cell[2][1]))
            fout.write('Masses\n\n  1   26.9815\n 2   47.8670\n\nAtoms\n\n')
            ii=1
            for i in range(int(self.atom_type_num[0])):
                fout.write('    '+str(int(ii))+'   1   '+str("%0.4f" % self.data[i][0])+'   '+str("%0.4f" % self.data[i][1])+'   '+str("%0.4f" % self.data[i][2])+'   '+str(self.data[i][3])+'   '+str(self.data[i][4])+'   '+str(self.data[i][5])+'\n')
                ii+=1
            for i in range(int(self.atom_type_num[0]),self.tot_num):
                fout.write('    '+str(int(ii))+'   2   '+str("%0.4f" % self.data[i][0])+'   '+str("%0.4f" % self.data[i][1])+'   '+str("%0.4f" % self.data[i][2])+'   '+str(self.data[i][3])+'   '+str(self.data[i][4])+'   '+str(self.data[i][5])+'\n')
                ii+=1
        return

    #Monte Carlo change two atoms
    def MC_ex(self):
        '''
        random change position of an atom with its neighbor
        '''
        a=rd.randint(0,self.tot_num-1)
        nb=[]
        for i in range(self.tot_num):
            if (calc_dist(self.data,a,i,self.cell) <= 3) and (i!=a):
                nb.append(i)
        #test use
        #print a,nb
        #for i in range(len(nb)):
        #    print calc_dist(self.data,a,nb[i],self.cell)
        bb=rd.randint(0,len(nb)-1)
        #print bb
        tmp_data=np.copy(self.data[nb[bb]])
        self.data[nb[bb]]=self.data[a]
        self.data[a]=tmp_data
        '''
        change between type 1 and type 2
        '''
        return

    def cal_EAM(self):
        os.system('mpirun -np 4 lmp_linux < lammps.in > md.out')
        with open('md.out','r') as fin1:
            for line in fin1:
                if 'Loop time of' in line:
                    stp_tot=int(line.split()[8])
                    break
        PotEng=[]
        with open('md.out','r') as fin1:
            for line in fin1:
                if 'Step PotEng' in line :
                    for i in range(stp_tot+1):
                        PotEng.append(float(fin1.next().strip().split()[1]))
                    break
        tmp=float(PotEng[-1])
        return  float(tmp)
    def cal_ANN(self):
        cwd=os.getcwd()
        cmdln='cp Ti.10t-10t.ann Al.10t-10t.ann '+cwd+str(self.count1)
        cwd=os.getcwd()
        with open('predict.EAM','w') as fout1:
            fout1.write('TYPES\n2\nTi\nAl\n\nNETWORKS\n Ti Ti.10t-10t.ann\n Al Al.10t-10t.ann\n\n!FORCES\n\nFILES\n1\n')
            fout1.write('now.xsf\n')
        os.system('mpirun -np 4 predict.mpi predict.EAM > predict.out')
        energy_A=[]
        atom_num=[]
        with open('predict.out','r') as fin1:
            for line in fin1:
                if 'Total energy' in line:
                    energy_A.append(float(line.split()[3]))
                if 'atoms' in line:
                    atom_num.append(float(line.split()[4]))
        return float(energy_A[-1])

    def __init__(self,filename,steps,method):
        self.filename=str(filename)
        self.tot_num=0
        self.cell=np.zeros((3,3))
        self.data=np.zeros((self.tot_num,6))
        self.tmp_data=np.zeros((self.tot_num,6))
        self.steps=int(steps)
        self.method=str(method)
        #self.atom_type_num comes from get_data() function
        self.get_data()

        self.lmp_new()
        self.energy=float(self.cal_EAM())

        #I/O files
        f0=open('all_steps.dat','w')
        f1=open('accepted_steps.dat','w')
        #main prog strucure 
        i0=0 #total steps to run
        i1=0 #total steps change structures
        print self.energy
        f0.write('%i    %12.6f\n'%(i0,self.energy))
        f1.write('%i    %12.6f  %i\n'%(i1,self.energy,i0))

        if self.method=='EAM':
            while (i0<=self.steps):
                self.tmp_data=np.copy(self.data)
                self.MC_ex()
                self.lmp_new()
                #if (i1%100 == 0):
                if (i1%1 == 0):
                    self.tmp_name=os.getcwd()+'/xsf/out_'+str(i1).zfill(6)+'.xsf'
                    self.xsf_new()
                energy_A=self.cal_EAM()
                if energy_A < self.energy:
                    self.energy=energy_A
                    i0+=1;i1+=1
                    f0.write('%i    %12.6f\n'%(i0,energy_A))
                    f1.write('%i    %12.6f  %i\n'%(i1,energy_A,i0))
                elif rd.random()<np.exp((self.energy-energy_A)/0.0259):
                    self.energy=energy_A
                    i0+=1;i1+=1           
                    f0.write('%i    %12.6f\n'%(i0,energy_A))
                    f1.write('%i    %12.6f  %i\n'%(i1,energy_A,i0))
                elif rd.random()>=np.exp((self.energy-energy_A)/0.0259):
                    i0+=1
                    f0.write('%i    %12.6f\n'%(i0,energy_A))
                    self.data=np.copy(self.tmp_data)

        elif  self.method=='ANN':
            while (i0<=self.steps):
                self.tmp_data=np.copy(self.data)
                self.MC_ex()
                self.tmp_name='now.xsf'
                self.xsf_new()
                #if (i1%100 == 0):
                if (i1%1 == 0):
                    self.tmp_name=os.getcwd()+'/xsf/out_'+str(i1).zfill(6)+'.xsf'
                    self.xsf_new()
                energy_A=self.cal_ANN()
                if energy_A < self.energy:
                    self.energy=energy_A
                    i0+=1;i1+=1
                    f0.write('%i    %12.6f\n'%(i0,energy_A))
                    f1.write('%i    %12.6f  %i\n'%(i1,energy_A,i0))
                elif rd.random()<np.exp((self.energy-energy_A)/0.0259):
                    self.energy=energy_A
                    i0+=1;i1+=1           
                    f0.write('%i    %12.6f\n'%(i0,energy_A))
                    f1.write('%i    %12.6f  %i\n'%(i1,energy_A,i0))
                elif rd.random()>=np.exp((self.energy-energy_A)/0.0259):
                    i0+=1
                    f0.write('%i    %12.6f\n'%(i0,energy_A))
                    self.data=np.copy(self.tmp_data)
        #main prog strucure end
        f0.close()
        f1.close()
        return

if __name__=="__main__":
    test=MonteCarlo('0_sf.xsf',1,'EAM')

