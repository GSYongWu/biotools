#!/usr/bin/env python

import multiprocessing as mp
import argparse
import sys
import fastq_infos
import numpy as np
import shlex
import subprocess
import modules


def get_fastqinfo(fastqlist):
    lis = [fastq_infos.fastqinfo(x) for x in fastqlist]
    dic = modules.combineFastqInfos(np.array(lis))
    return dic



