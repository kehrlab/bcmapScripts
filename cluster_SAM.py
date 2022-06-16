# cluster alignemnt data in SAM format into barcode mappings
# the input file has to be sorted by barcode

import sys
from itertools import islice

class posSet:
    def __init__(self,pos):
        self.positions=[pos]

    def insert(self,pos):
        i=0
        for po in self.positions:
            if pos.smallerThan(po):
                self.positions.insert(i,pos)
                return
            i+=1
        self.positions.append(pos)
        return

class position:
    def __init__(self,chr,pos1,pos2):
        self.chr=chr
        if pos1<pos2:
            self.pos1=pos1
            self.pos2=pos2
        else:
            self.pos1=pos2
            self.pos2=pos1

    def smallerThan(self,pos):
        if self.chr<pos.chr:
            return 1
        elif self.chr>pos.chr:
            return 0
        if int(self.pos1)<int(pos.pos1):
            return 1
        elif int(self.pos1)>int(pos.pos1):
            return 0
        if int(self.pos2)<int(pos.pos2):
            return 1
        else:
            return 0

class mapping:
    def __init__(self,pos):
        self.chr=pos.chr
        self.start=pos.pos1
        self.end=pos.pos2
        self.reads=1

def getposition(line):
    splitline=line.split('\t')
    pos=position(splitline[2],splitline[3],splitline[7])
    # print(pos.chr, ' ', pos.pos1, ' ', pos.pos2)
    return pos

def getreadbarcode(line):
    barcode=line.split()[1][5:]
    # print("barcode: ",barcode)
    return barcode

if len(sys.argv)!=4:
    print('usage: cluster_SAM.py sim_readfile1.fq alignment.sam outputfile.bed')
    print('')
    print('Cluster alignemnt data in SAM format into barcode mappings')
    print('The input file has to be sorted by barcode')
    exit()

# define the number of lines to read
old_barcode='NNN'
number_of_lines=1000
max_gap=20000
max_len=300000
min_len=10000
barcode_count=0
# min_reads_per_mapping=int(sys.argv[4]) # 10
map=mapping(position('nochr','0','0'))
outputfile=open(sys.argv[3],'w')
samfile=open(sys.argv[2],'r')
mappings=''
positions=posSet(position('nochr','0','0'))
with open(sys.argv[1],'r') as input_file:
    while 1:
        lines_cache=list(islice(input_file, number_of_lines))
        nol=int(number_of_lines/2)
        sam_cache=list(islice(samfile,nol))
        if not lines_cache:
            if abs(int(map.end)-int(map.start))>min_len:
                mappings+=map.chr+'\t'+map.start+'\t'+map.end+'\t'+old_barcode+'\t'+str(map.reads)+'\n'
            outputfile.write(mappings)
            print('barcodes: ',barcode_count )
            print('i managed to break!')
            break
        for i in range(0,len(lines_cache),4):
            j=i/4
            barcode=getreadbarcode(lines_cache[i])
            # for j in range(int(i/2),int(i/2+2)):
                # print("i: ",i," j: ",j)
            pos=getposition(sam_cache[int(i/2)])
            if barcode==old_barcode:
                # append posSet
                # print("Pos: [",pos.chr,", ",pos.pos1,", ",pos.pos2)
                positions.insert(pos)
            else:
                barcode_count+=1
                #evaluate posset
                # print("posset: ")
                # if old_barcode=="AAACACCAGCATGTTC":
                #     for posi in positions.positions:
                #         print("[",posi.chr,", ",posi.pos1,", ",posi.pos2)

                for posi in positions.positions:
                    if posi.chr==map.chr and abs(int(posi.pos1)-int(map.end))<max_gap and abs(int(posi.pos2)-int(map.start))<max_len:
                        map.end=posi.pos2
                        map.reads+=1
                    else:
                        if abs(int(map.end)-int(map.start))>min_len:
                            mappings+=map.chr+'\t'+map.start+'\t'+map.end+'\t'+old_barcode+'\t'+str(map.reads)+'\n'
                            if len(mappings)>10000:
                                outputfile.write(mappings)
                                mappings=''
                        map=mapping(posi)
                #report last map
                if abs(int(map.end)-int(map.start))>min_len:
                    mappings+=map.chr+'\t'+map.start+'\t'+map.end+'\t'+old_barcode+'\t'+str(map.reads)+'\n'
                    if len(mappings)>10000:
                        outputfile.write(mappings)
                        mappings=''
                #create new posset
                positions=posSet(pos)
                old_barcode=barcode
                map=mapping(position('nochr','0','0'))
