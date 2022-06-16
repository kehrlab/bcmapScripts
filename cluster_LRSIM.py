# Read data simulated with LRSIM has to be barcode trimmed and sorted with bcctools.
# Use this script to cluster the corrected fastq files into a barcode mapping bed file serving as a truth set.
# remove unidentified barcodes marked by a * (manually)
# Compare barcode mappings to the truth set using the script: "compare_barcode_mappings.py".

import sys
from itertools import islice

class position:
    def __init__(self,chr,pos1,pos2):
        self.chr=chr
        if pos1<pos2:
            self.pos1=pos1
            self.pos2=pos2
        else:
            self.pos1=pos2
            self.pos2=pos1

class mapping:
    def __init__(self,pos):
        self.chr=pos.chr
        self.start=pos.pos1
        self.end=pos.pos2

def getposition(line):
    splitline=line.split('_')[0:3]
    pos=position(splitline[0][1:],splitline[1],pos2=splitline[2])
    # print(pos.chr, ' ', pos.pos1, ' ', pos.pos2)
    return pos

def getreadbarcode(line):
    barcode=line.split()[1][5:]
    # print("barcode: ",barcode)
    return barcode

if len(sys.argv)!=3:
    print('usage: cluster_LRSIM.py sim_readfile1.fq outputfile.bed')
    print('')
    print('Read data simulated with LRSIM has to be barcode trimmed and sorted with bcctools.')
    print('Then use this script to cluster the corrected fastq files into a barcode mapping bed file serving as a truth set.')
    print('Remove unidentified barcodes marked by a * (manually)')
    print('Compare barcode mappings to the truth set using the script: "compare_barcode_mappings.py".')
    exit()

# define the number of lines to read
old_barcode='*'
number_of_lines=1000
max_gap=20000
max_len=300000
min_len=10000
map=mapping(position('nochr','0','0'))
outputfile=open(sys.argv[2],'w')
mappings=''
with open(sys.argv[1],'r') as input_file:
    while 1:
        lines_cache=list(islice(input_file, number_of_lines))
        if not lines_cache:
            if abs(int(map.end)-int(map.start))>min_len:
                mappings+=map.chr+'\t'+map.start+'\t'+map.end+'\t'+old_barcode+'\n'
            outputfile.write(mappings)
            print('i managed to break!')
            break
        for i in range(0,len(lines_cache),4):
            barcode=getreadbarcode(lines_cache[i])
            pos=getposition(lines_cache[i])
            if barcode==old_barcode:
                if pos.chr==map.chr and abs(int(pos.pos1)-int(map.end))<max_gap and abs(int(pos.pos2)-int(map.start))<max_len:
                    map.end=pos.pos2
                else:
                    if abs(int(map.end)-int(map.start))>min_len:
                        mappings+=map.chr+'\t'+map.start+'\t'+map.end+'\t'+barcode+'\n'
                        if len(mappings)>10000:
                            outputfile.write(mappings)
                            mappings=''
-                    map=mapping(pos)
            else:
                if abs(int(map.end)-int(map.start))>min_len:
                    mappings+=map.chr+'\t'+map.start+'\t'+map.end+'\t'+old_barcode+'\n'
                    if len(mappings)>10000:
                        outputfile.write(mappings)
                        mappings=''
                map=mapping(pos)
                old_barcode=barcode
