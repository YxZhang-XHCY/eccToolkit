#coding=utf-8
import sys
aList=[]
fa_file = sys.argv[1]
with open(fa_file,'r') as f:
    for line in f:
        line = line.strip()
        line = line.upper()
        if not line.startswith(">"):
            baseA = line.count("A")
            baseT = line.count("T")
            baseC = line.count("C")
            baseG = line.count("G")
            aList.extend([baseA, baseT, baseC, baseG])
            # print(aList)
    print("effective_genome_size =", sum(aList))
