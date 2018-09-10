#!/usr/bin/env python
import os
import subprocess
import threading

from pybedtools import BedTool


def bamdrs_run(bamdrspath,
               sortbam,
               rmdupbam,
               bed,
               outdir,
               flank=100,
               mapQ=20,
               uncover=20):
    rawdir = outdir + "/raw"
    flankdir = outdir + "/flank"
    rmdupdir = outdir + "/rmdup"
    flankbed = flankdir + "/flank.bed"
    makedir(rawdir, flankdir, rmdupdir)
    slopbed(bed, flank, flankbed)

    threads = []
    raw_process = threading.Thread(
        target=bamdrs_subprocess,
        args=(
            bamdrspath,
            sortbam,
            bed,
            rawdir,
            mapQ,
            uncover,
        ))
    flank_process = threading.Thread(
        target=bamdrs_subprocess,
        args=(
            bamdrspath,
            sortbam,
            flankbed,
            flankdir,
            mapQ,
            uncover,
        ))
    rmdup_process = threading.Thread(
        target=bamdrs_subprocess,
        args=(
            bamdrspath,
            rmdupbam,
            bed,
            rmdupdir,
            mapQ,
            uncover,
        ))
    threads = [raw_process, flank_process, rmdup_process]
    for t in threads:
        t.setDaemon(True)
        t.start()
    t.join()


def bamdrs_subprocess(bamdrspath, bam, bed, outdir, mapQ=20, uncover=20):
    """Generate bamdrs command list
    
    Arguments:
        bamdrspath {[string]} -- [path to bamdrs]
        bam {[string]} -- [path to bam file]
        bed {[string]} -- [path to bed file]
        outdir {[string]} -- [path to outdir]
    
    Keyword Arguments:
        mapQ {int} -- [description] (default: {20})
        uncover {int} -- [description] (default: {20})
    
    Returns:
        [list] -- [command line list]
    """

    cmdlis = [
        bamdrspath, '-q', mapQ, '--uncover', uncover, '-p', bed, '-o', outdir,
        bam
    ]
    subprocess.call(cmdlis)


def makedir(*dirs):
    for path in dirs:
        if not os.path.exists(path):
            os.makedirs(path)


def slopbed(bedfile, flank, outfile):
    """Extend bed file.
    
    Arguments:
        bedfile {[string]} -- [raw bed file]
        flank {[int]} -- [extend base pairs in each direction]
        outfile {[string]} -- [extend bed file save as]
    """

    rawbed = BedTool(bedfile)
    addbed = rawbed.slop(b=flank, genome="hg19").sort().merge()
    addbed.saveas(outfile)


def coverage2dict(coverage_report):
    """bamdst coverage.report file to dict.
    
    Arguments:
        coverage_report {[string]} -- [bamdst coverage.report file.]
    
    Returns:
        [dict] -- [dict of coverage.report infos]
    """

    dic = {}
    with open(coverage_report) as fin:
        for i in fin:
            if i.startswith("#"):
                continue
            i = i.strip()
            t = i.split("\t")
            dic[t[0]] = t[1]
    return dic


def nxcoverage(depth_distribution, ndepth=500):
    """Calculate coverage when depth >=N and median depth.
    
    Arguments:
        depth_distribution {string} -- [bamdst output depth_distribution.plot file]
    
    Keyword Arguments:
        ndepth {int} -- [depth for calculate NXcoverage] (default: {500})

    Returns:
        [tuple] -- [NXdepthCoverage and median depth]
    """

    median_depth = 0
    coverage = 0
    #depth_lis = []
    coverage_lis = []
    with open(depth_distribution) as fin:
        for i in fin:
            i = i.rstrip()
            depth, base, freq, rest_base, sumfreq = i.split("\t")
            coverage_lis.append(float(sumfreq))
    coverage = coverage_lis[ndepth - 1]
    for i in range(len(coverage_lis) + 1):
        if coverage_lis[i] > 0.5 and coverage_lis[i + 1] <= 0.5:
            median_depth = i + 1
            break
    return coverage, median_depth


def insert_size(insertsize_plot):
    """Calculate median insert size
    
    Arguments:
        insertsize_plot {[string]} -- [bamdst insertsize.plot output file ]
    
    Returns:
        [int] -- [median insert size]
    """

    median_insert = 0
    freq_lis = []
    with open(insertsize_plot) as fin:
        for i in fin:
            i = i.rstrip()
            t = i.split("\t")
            freq_lis.append(float(t[-1]))
    for i in range(len(freq_lis) + 1):
        if freq_lis[i] > 0.5 and freq_lis[i + 1] <= 0.5:
            median_insert = i + 1
            break
    return median_insert


if __name__ == '__main__':

    import sys
    if len(sys.argv) < 3:
        print "py bamdrs sortbam rmdupbam bed outdir"
        exit(1)
    bamdrs_run(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
