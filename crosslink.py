#!/usr/bin/env python
import sys
import os
import re
from collections import defaultdict
"""
1. mask read mismatch >= N1
2. mask read with clip and mismatch
3. mask read tail
"""
if len(sys.argv) < 2:
    print "py in.bam /path/to/samtools hq_bam.log"
    exit(1)
inbam = sys.argv[1]
samtools = sys.argv[2]
log_file = open(sys.argv[3], 'w')

mismatch = 5
total_read = 0
pass_reads = 0
filtered_reads = 0
for i in os.popen('{samtools} view -h {inbam}'.format(
        samtools=samtools, inbam=inbam)):
    #print i,
    if i.startswith('@'):
        print i,
    else:
        total_read += 1
        t = i.split('\t')
        mapq = int(t[4])
        #readid = t[0]
        #chrom = t[2]
        #start = int(t[3])
        #cigar = t[5]
        #cigarlist=re.findall("[1-9][0-9]*[MIDSHNP=X]",cigar)
        #cigardict=defaultdict(lambda:0)
        #for n in cigarlist:
        #	cigardict[n[-1]]+=int(n[:-1])
        #reflen=cigardict['M']+cigardict['D']
        #readlen=cigardict['M']+cigardict['I']+cigardict['S']+cigardict['H']
        #end = start + readlen

        flag = int(t[1])
        if flag & 256 == 256:
            primary = False
        elif flag & 512 == 512:
            primary = False
        elif flag & 1024 == 1024:
            primary = False
        elif flag & 2048 == 2048:
            primary = False
        else:
            primary = True

        cigar = t[5]
        clip = 0
        if 'S' in cigar:
            clip = int(max(re.findall(r'(\d+)S', cigar)))
        nins = cigar.count('I')
        nsnp = 0
        ndel = 0
        mdinfos = re.findall(r'\tMD:Z:(\S+)\t', i)
        if mdinfos:
            mdinfo = mdinfos[0]
            muts = re.findall(r'[ATCGNatcgn]|\^[ATCGNatcgn]+', mdinfo)
            for mut in muts:
                if len(mut) == 1:
                    nsnp += 1
                else:
                    ndel += 1
        N = nins + nsnp + ndel
        if clip:
            N += 1
        #Mismatch check
        if primary and N < mismatch and mapq > 30:
            pass_reads += 1
            print i,
        else:
            filtered_reads += 1

pass_reads_pct = pass_reads * 100.0 / total_read
filtered_reads_pct = filtered_reads * 100.0 / total_read
sid = inbam.split('/')[-1].split('.')[0]
log_file.writelines('#ID,PassReadsPct,FailReadsPct\n{},{},{}\n'.format(
    sid, round(pass_reads_pct, 2), round(filtered_reads_pct, 2)))
log_file.close()


""" /usr/bin/python /GPFS01/databases/GSCAP2018/Scripts/misc/CleanBam.py 
/GPFS01/AnalysisTemp/pipeline2018_test/INT/Results/B999/TMP/B999.sorted.rmdup.bam 
 /GPFS01/databases/GSCAP2018/softwares/samtools-1.2/samtools 
 /GPFS01/AnalysisTemp/pipeline2018_test/INT/Results/B999/QCSummary/B999.clean.bam.csv 
 >/dev/null 2>&1

"""
