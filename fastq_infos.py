#!/usr/bin/env python
from __future__ import print_function

import os
import Queue
import subprocess
import sys
import threading


def basecount(seqstring):
    GC = seqstring.count("G") + seqstring.count("C")
    N = seqstring.count("N")
    return len(seqstring), GC, N


def qualitycount(qualstring):
    #chr(20+33)='5' chr(30+33)='?'
    q20 = 0
    q30 = 0
    for n in qualstring:
        if n >= '?':
            q20 += 1
            q30 += 1
        elif n >= '5':
            q20 += 1
    return q20, q30


def fastqreader(fastqinpath):
    if fastqinpath == '-':
        filein = sys.stdin
    else:
        if fastqinpath[-9:] == ".fastq.gz" or fastqinpath[-6:] == ".fq.gz":
            gzfastqin = subprocess.Popen(
                ('unpigz', '-c', fastqinpath), stdout=subprocess.PIPE)
            filein = gzfastqin.stdout
        else:
            filein = open(fastqinpath)
        #elif fastqinpath[-6:]==".fastq" or fastqinpath[-3:]==".fq":
        #    gzfastqin=subprocess.Popen(('cat',fastqinpath),stdout=subprocess.PIPE)
        #    filein=gzfastqin.stdout
        #else:
        #    raw_input("Input Fastq file format is ERROR! %s" %fastqinpath)
        #    sys.exit()
    count = 0
    for line in filein:
        if line != '':
            count += 1
            #get the second line of a read (seq line)
            seqlen, gccount, ncount = basecount(filein.next().strip())
            filein.next()
            #get the fourth line of read (qual line)
            q20, q30 = qualitycount(filein.next().strip())
            yield seqlen, gccount, ncount, q20, q30
        else:
            break
    filein.close()


def fastqinfo(fastqfile):
    """Calculate fatsq infos
    
    Arguments:
        fastqfile {[str]} -- [path to fastq file]
    """

    filename = os.path.basename(fastqfile)
    readnum = 0
    maxlen = 0
    basenum = 0
    GCcount = 0
    Ncount = 0
    Q20 = 0
    Q30 = 0
    fastqin = fastqreader(fastqfile)
    for seqlen, gccount, ncount, q20, q30 in fastqin:
        readnum += 1
        maxlen = max(maxlen, seqlen)
        basenum += seqlen
        GCcount += gccount
        Ncount += ncount
        Q20 += q20
        Q30 += q30


#    print('#SAMPLE,READS,LEGTH,BASES,GC,Q20,Q30,PPM')

    return [
        filename, readnum, maxlen,
        max(round(basenum / 1000 / 1000., 2), 0.01),
        round(GCcount * 100. / basenum, 2),
        round(Q20 * 100. / basenum, 2),
        round(Q30 * 100. / basenum, 2),
        round(Ncount * 1000. * 1000. / basenum, 2)
    ]

"""
def combineR1R2(lisinfo1, lisinfo2):
    dic = {}
    filename = lisinfo1[0]
    dic["SAMPLE"] = filename.split("_")[0]
    dic["PF_READS"] = lisinfo1[1]
    dic["LENGTH"] = lisinfo1[2]
    dic["PF_BASES"] = lisinfo1[3] + lisinfo2[2]
    dic["GC(%)"] = round((lisinfo1[4] + lisinfo2[4]) / 2, 2)
    dic["Q20(%)"] = round((lisinfo1[5] + lisinfo2[5]) / 2, 2)
    dic["Q30(%)"] = round((lisinfo1[6] + lisinfo2[6]) / 2, 2)
    return dic


def thread_process(fastqR1, fastqR2):
    q = Queue.Queue()

    def dealforqueue(fastqfile):
        q.put(fastqinfo(fastqfile))

    result = []
    t1 = threading.Thread(
        target=dealforqueue, name="fastqR1", args=(fastqR1, ))
    t2 = threading.Thread(
        target=dealforqueue, name="fastqR2", args=(fastqR2, ))
    t1.start()
    t2.start()
    t1.join()
    t2.join()
    while not q.empty():
        result.append(q.get())
    dic = combineR1R2(result[0], result[1])
    return dic

"""
def main():
    infolis = []
    if len(sys.argv) == 2:
        assert os.path.exists(
            sys.argv[1]
        ), "File %s does not exist.\nUsage:\npython script.py filepath\npython - sampleid" % sys.argv[1]
        infolis = fastqinfo(sys.argv[1])
    else:
        print("Usage:\npython script.py filepath\npython script.py R1.fastq.gz R2.fastq.gz")
    if infolis:
        #print('#SAMPLE,READS,LEGTH,BASES,GC,Q20,Q30,PPM')
        print(','.join(map(str, infolis)))


if __name__ == "__main__":
    main()
