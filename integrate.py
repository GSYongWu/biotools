#!/usr/bin/env python

import sys
import fastq_infos
import numpy as np
from bamdst_parse  import bamdst_run
from modules import *
import multiprocessing as mp
from options import *

def get_fastqinfo(input_fastqs):
    list_fastq, list_lane = raw_fastq(input_fastqs)
    lis = [fastq_infos.fastqinfo(x) for x in list_fastq]
    dic = combineFastqInfos(np.array(lis))
    return dic


def program():

    pass

def run_all(config, inputfastqs, sortbam, rmdupbam, bedfile ):
    cfg = readconfig(config)
    pool = mp.Pool(processes=6)
    results = []
    results.append(
        pool.apply_async(get_fastqinfo, (inputfastqs, ))
    )
    results.append(
        pool.apply_async(bamdst_run, (cfg.get("software", "bamdst"),
        sortbam, rmdupbam, bedfile, cfg.get("params", "genome")+".fai",cfg.get("software","bedtools") )),
    )
    output = [p.get() for p in results]
    dic = {}
    for d in output:
        dic.update(d)
    return dic


def main():
    args = parseargs()
    run_all(args.config, args.fastqinput, args.sortedbam, args.rmdupbam, args.bed)



if __name__ == '__main__':
    main()