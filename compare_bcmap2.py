import sys
import numpy as np
import math

class mapping:
    def __init__(self, line):
        split=line.split()
        # print(split)
        self.chr=split[0]
        self.start=int(split[1])
        self.end=int(split[2])
        self.barcode=split[3]
        self.score=0

    def show(self):
        print(self.chr,'\t',self.start,'\t',self.end,'\t',self.barcode)


class bc_mapping:
    def __init__(self, line):
        split=line.split()
        # print(split)
        self.chr=split[0]
        self.start=int(split[1])
        self.end=int(split[2])
        self.barcode=split[3]
        self.score=float(split[4])
        # self.score=int(split[4])/(self.end-self.start)

    def show(self):
        print(self.chr,'\t',self.start,'\t',self.end,'\t',self.barcode, '\t',self.score)

class result:
    def __init__(self):
        # self.score=[0.001,0.002,0.004,0.006,0.008,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1]
        # self.score=[0.01,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3]
        # self.score=[0.01,0.05,0.1,0.125,0.15,0.175,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,1,1.5,2.0,2.5]
        self.score=[0.5,1,1.5,2,2.5,3,3.5,4,6,8,10,12,14,16,18,20,22,24,26,28,30,40,50,60,70,80,90,100,120,140,160,180,200]
        # self.score=[30,40,50,100,200,300,400,500,600,700,800,900,1000,1200,1400,1600,1800,2000]
        # self.score=[50,500,1000,2000,3000,4000,5000,10000,15000,20000,25000,30000,31000,32000,33000,34000,35000,40000,50000,100000]
        self.TP=np.zeros(len(self.score),dtype=int)
        self.FP=np.zeros(len(self.score),dtype=int)
        self.FN=np.zeros(len(self.score),dtype=int)
        self.DEV=np.zeros(len(self.score),dtype=int)

def getreadbarcode(line):
    barcode=line.split()[1][5:]
    # print("barcode: ",barcode)
    return barcode

def getposition(line):
    splitline=line.split('_')[0:3]
    pos=position(splitline[0][1:],splitline[1],pos2=splitline[2])
    # print(pos.chr, ' ', pos.pos1, ' ', pos.pos2)
    return pos

# version 2
# def match(bcmap,lrsim):
#     if bcmap.barcode!=lrsim.barcode:
#         print('Hey, this is weird!')
#     if bcmap.chr!=lrsim.chr:
#         return 0
#     if bcmap.start>lrsim.start and bcmap.start<lrsim.end: # lrsim contains start of bcmap
#         return 1
#     if bcmap.end>lrsim.start and bcmap.end<lrsim.end: # lrsim contains end of bcmap
#         return 1
#     if bcmap.start<lrsim.start and bcmap.end>lrsim.end: # bcmap contains lrsim
#         return 1
#     return 0

# version1
def match(bcmap,lrsim):
    max_dev=20000
    if bcmap.barcode!=lrsim.barcode:
        print('Hey, this is weird!')
    if bcmap.chr!=lrsim.chr:
        return 0
    if abs(bcmap.start-lrsim.start)>max_dev:
        return 0
    if abs(bcmap.end-lrsim.end)>max_dev:
        return 0
    if bcmap.start<lrsim.start and bcmap.end<lrsim.start:
        return 0
    if bcmap.end>lrsim.end and bcmap.start>lrsim.end:
        return 0
    return 1

def compare(bcmap_mappings,lrsim_mappings,res):
    for bcmap in bcmap_mappings:
        for lrsim in lrsim_mappings:
            if match(bcmap,lrsim):
                for i in range(len(res.score)):
                    if bcmap.score>res.score[i]:
                        if res.score[i]>lrsim.score:
                            lrsim.score=res.score[i]
                            res.DEV[i]+=abs(bcmap.end-lrsim.end)+abs(bcmap.start-lrsim.start)
                            res.TP[i]+=1
                        else:
                            res.FP[i]+=1
                break
        else:
            for i in range(len(res.score)):
                if bcmap.score>res.score[i]:
                    res.FP[i]+=1
                    # if i==17:
                        # bcmap.show()
    for lrsim in lrsim_mappings:
        for i in range(len(res.score)):
            if res.score[i] > lrsim.score:
                res.FN[i]+=1

    # for lrsim in lrsim_mappings:
    #     for bcmap in bcmap_mappings:
    #         if match(bcmap,lrsim):
    #             break
    #     else:
    #         # lrsim.show()
    #         for i in range(len(res.score)):
    #             if bcmap.score>res.score[i]:
    #                 res.FN[i]+=1
    return res

if len(sys.argv)!=3:
    print('usage: compare_bcmap2.py res_sim.bed res_bcmap.bed')
    exit()

sim_res=open(sys.argv[1],'r')
bcmap_res=open(sys.argv[2],'r')
bcmap_mapping=bc_mapping(bcmap_res.readline())
bcmap_mappings=[bcmap_mapping]
lrsim_mapping=mapping(sim_res.readline())
lrsim_mappings=[lrsim_mapping]
res=result()

if bcmap_mapping.barcode!=lrsim_mapping.barcode:
    print('ERROR! HELP ME!')
else:
    barcode=bcmap_mapping.barcode


for line in bcmap_res:
    bcmap_mapping=bc_mapping(line)
    if bcmap_mapping.barcode==barcode:
        bcmap_mappings.append(bcmap_mapping)
    else:
        lrsline=sim_res.readline()
        if lrsline=='':
            print('TP: ',res.TP,' FP: ',res.FP,' FN: ',res.FN, 'DEV: ',res.DEV/(res.TP*2))
            exit()
        lrsim_mapping=mapping(lrsline)
        while lrsim_mapping.barcode==barcode:
            lrsim_mappings.append(lrsim_mapping)
            lrsline=sim_res.readline()
            if lrsline=='':
                res=compare(bcmap_mappings,lrsim_mappings,res)
                print('TP: ',res.TP,' FP: ',res.FP,' FN: ',res.FN)
                exit()
            lrsim_mapping=mapping(lrsline)
        # comparison
        res=compare(bcmap_mappings,lrsim_mappings,res)
        # reinitialization
        while lrsim_mapping.barcode!=bcmap_mapping.barcode:
            while lrsim_mapping.barcode<bcmap_mapping.barcode:
                lrsline=sim_res.readline()
                for i in range(len(res.score)):
                    res.FN[i]+=1
                if lrsline=='':
                    print('TP: ',res.TP,' FP: ',res.FP,' FN: ',res.FN)
                    exit()
                lrsim_mapping=mapping(lrsline)
            while lrsim_mapping.barcode>bcmap_mapping.barcode:
                for i in range(len(res.score)):
                    res.FP[i]+=1
                # print('Hi!')
                line=next(bcmap_res)
                if line!='':
                    bcmap_mapping=bc_mapping(line)
                else:
                    print('Help, i need somebody, Help!')
            break
                # print(lrsim_mapping.barcode)
                # print(bcmap_mapping.barcode)
        barcode=bcmap_mapping.barcode
        bcmap_mappings=[bcmap_mapping]
        lrsim_mappings=[lrsim_mapping]

# print('Score: ',res.score,' \nTP: ',res.TP,' \nFP: ',res.FP,' \nFN: ',res.FN)
Precision=[round(res.TP[i]/(res.TP[i]+res.FP[i])*100,2) for i in range(len(res.score)) if res.TP[i]>0]
Recall=[round(res.TP[i]/(res.TP[i]+res.FN[i])*100,2) for i in range(len(res.score)) if res.TP[i]>0]
Deviation=[round(res.DEV[i]/(res.TP*2)[i]) for i in range(len(res.score)) if res.TP[i]>0]
print('Precision: ', Precision, ' \nRecall: ', Recall, '\nDEV: ', Deviation)
