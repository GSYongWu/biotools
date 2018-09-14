#!/usr/bin/env python
import argparse
import gzip
import os
import subprocess
import sys
import threading
from time import ctime

import numpy

from pybedtools import BedTool


def makedir(*dirs):
    for path in dirs:
        if not os.path.exists(path):
            os.makedirs(path)


def str2pctstr(floatstr):
    return str(round(float(floatstr) * 100, 2)) + "%"


# def slopbed(bedfile, flank, outfile):
#     """Extend bed file.

#     Arguments:
#         bedfile {[str]} -- [raw bed file]
#         flank {[int]} -- [extend base pairs in each direction]
#         outfile {[str]} -- [extend bed file save as]
#     """

#     rawbed = BedTool(bedfile)
#     addbed = rawbed.slop(b=flank, genome="hg19").sort().merge()
#     addbed.saveas(outfile)


def slopbed(bedfile, flank, outfile, genome, bedtools):
    cmd = "{bedtools} slop -i {bed} -b {flank} -g {genome} |"\
    "bedtools sort -i - | bedtools merge -i - > {outfile}".format(bedtools=bedtools,
     bed=bedfile, flank=flank, genome=genome, outfile=outfile)
    subprocess.Popen(cmd, shell=True)


def bamdst_subprocess(bamdstpath, bam, bed, outdir, mapQ=20, uncover=20):
    """Run bamdst command in subprocess

    Arguments:
        bamdstpath {[str]} -- [path to bamdst]
        bam {[str]} -- [path to bam file]
        bed {[str]} -- [path to bed file]
        outdir {[str]} -- [path to outdir]

    Keyword Arguments:
        mapQ {int} -- [map qual] (default: {20})
        uncover {int} -- [description] (default: {20})
    """

    cmdlis = [
        bamdstpath, '-q',
        str(mapQ), '--uncover',
        str(uncover), '-p', bed, '-o', outdir, bam
    ]
    subprocess.call(cmdlis)


def coverage2dict(coverage_report):
    """bamdst coverage.report file to dict.

    Arguments:
        coverage_report {[str]} -- [bamdst coverage.report file.]

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


def nxcoverage(depth_distribution, *ndepth):
    """Calculate coverage when depth >=N and median depth.

    Arguments:
        depth_distribution {str} -- [bamdst output depth_distribution.plot file]

    Returns:
        [tuple] -- [median depth and NXdepthCoverage list]
    """

    median_depth = 0
    coverage = []
    coverage_lis = []
    with open(depth_distribution) as fin:
        for i in fin:
            i = i.rstrip()
            depth, base, freq, rest_base, sumfreq = i.split("\t")
            coverage_lis.append(float(sumfreq))
    for nd in ndepth:
        coverage.append(str2pctstr(coverage_lis[nd - 1]))
    for i in range(len(coverage_lis) + 1):
        if coverage_lis[i] > 0.5 and coverage_lis[i + 1] <= 0.5:
            median_depth = i + 1
            break
    return median_depth, coverage


def insert_size(insertsize_plot):
    """Calculate median insert size

    Arguments:
        insertsize_plot {[str]} -- [bamdst insertsize.plot output file ]

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


def calCV(region_tsv_gz, stat="mean"):
    """Calculate CV score
    
    Arguments:
        region_tsv_gz {[str]} -- [bamdst output region.tsv.gz]
    
    Keyword Arguments:
        stat {str} -- [region depth stat method[median or mean]] (default: {"mean"})
    
    Returns:
        [float] -- [cv score]
    """

    mean_depth_lis = []
    median_depth_lis = []
    with gzip.open(region_tsv_gz) as fin:
        for i in fin:
            if i.startswith("#"):
                continue
            i = i.rstrip()
            t = i.split()
            mean_depth_lis.append(float(t[3]))
            median_depth_lis.append(float(t[4]))
    depthlis = []
    if stat == "mean":
        depthlis = mean_depth_lis
    elif stat == "median":
        depthlis = median_depth_lis
    else:
        sys.stderr.write("stat [mean or median]")
    first_quartile = numpy.percentile(depthlis, 25)
    third_quartile = numpy.percentile(depthlis, 75)
    iqr = third_quartile - first_quartile
    median_depth = numpy.median(depthlis)
    if median_depth == 0:
        iqr_median_cv = 999
    else:
        iqr_median_cv = iqr * 1.0 / median_depth
    return round(iqr_median_cv, 2)


def bamdst_integrate(sampleid, coverage_sort, depth_distribution_sort,
                     insertsize_sort, region_sort, coverage_rmdup,
                     depth_distribution_rmdup, coverage_flank):
    """Integrate bamdst output(sort bam raw bed / sort bam flank bed / rmdup bam raw bed)

    Arguments:
        sampleid {[str]} -- [sample id]
        coverage_sort {[str]} -- [path to sort coverage.report]
        depth_distribution_sort {[str]} -- [path to sort depth_distribution.plot]
        insertsize_sort {[str]} -- [path to sort insertsize.plot]
        coverage_rmdup {[str]} -- [path to rmdup coverage.report]
        depth_distribution_rmdup {[str]} -- [path to rmdup depth_distribution.plot]
        coverage_flank {[str]} -- [path to flank coverage.report]

    Returns:
        [dict] -- [sample infos]
    """

    dic = {}
    dic_sort = coverage2dict(coverage_sort)
    dic_rmdup = coverage2dict(coverage_rmdup)
    dic_flank = coverage2dict(coverage_flank)
    dic["#SAMPLE"] = sampleid
    dic["CLEAN_READS"] = int(
        dic_sort["[Total] Raw Reads (All reads)"]) / 2  #int(read pair)
    dic["CLEAN_BASES"] = dic_sort["[Total] Raw Data(Mb)"]  #MB, string
    dic["MAPQ20_READS"] = dic_sort["[Total] MapQuality above cutoff reads"]
    dic["MAPQ20(%)"] = dic_sort[
        "[Total] Fraction of MapQ reads in mapped reads"]
    dic["TARGET_READS"] = dic_sort["[Target] Target Reads"]
    dic["MAPPED_READS"] = dic_sort["[Total] Mapped Reads"]
    dic["ON_TARGET_CORE(%)"] = dic_sort[
        "[Target] Fraction of Target Data in all data"]
    dic["ON_TARGET_READS_CORE(%)"] = dic_sort[
        "[Target] Fraction of Target Reads in all reads"]
    dic["RATIO_OF_MAPPED(%)"] = dic_sort["[Total] Fraction of Mapped Reads"]
    dic["MEAN_DEPTH"] = float(dic_sort["[Target] Average depth"])
    dic["1X_COVERAGE(%)"] = dic_sort["[Target] Coverage (>0x)"]
    dic["100X_COVERAGE(%)"] = dic_sort["[Target] Coverage (>=100x)"]
    dic["INSERT_SIZE"] = insert_size(insertsize_sort)

    depth10pct, depth20pct, depth50pct = int(dic["MEAN_DEPTH"]) / 10, int(
        dic["MEAN_DEPTH"]) / 5, int(dic["MEAN_DEPTH"]) / 2

    dic["MEDIAN_DEPTH"], [
        dic["20X_COVERAGE(%)"], dic["50X_COVERAGE(%)"],
        dic["200X_COVERAGE(%)"], dic["500X_COVERAGE(%)"],
        dic["10%_COVERAGE(%)"], dic["20%_COVERAGE(%)"], dic["50%_COVERAGE(%)"]
    ] = nxcoverage(depth_distribution_sort, 20, 50, 200, 500, depth10pct,
                   depth20pct, depth50pct)

    dic["MEAN_DEPTH_DEDUP"] = float(dic_rmdup["[Target] Average depth"])
    dic["1X_COVERAGE_DEDUP(%)"] = dic_rmdup["[Target] Coverage (>0x)"]
    dic["100X_COVERAGE_DEDUP(%)"] = dic_rmdup["[Target] Coverage (>=100x)"]
    depthdup10pct, depthdup20pct, depthdup50pct = int(
        dic["MEAN_DEPTH_DEDUP"]) / 10, int(dic["MEAN_DEPTH_DEDUP"]) / 5, int(
            dic["MEAN_DEPTH_DEDUP"]) / 2
    dic["MEDIAN_DEPTH_DEDUP"], [
        dic["20X_COVERAGE_DEDUP(%)"], dic["50X_COVERAGE_DEDUP(%)"],
        dic["200X_COVERAGE_DEDUP(%)"], dic["500X_COVERAGE_DEDUP(%)"],
        dic["10%MEAN_COVERAGE_DEDUP(%)"], dic["20%MEAN_COVERAGE_DEDUP(%)"],
        dic["50%MEAN_COVERAGE_DEDUP(%)"]
    ] = nxcoverage(depth_distribution_rmdup, 20, 50, 200, 500, depthdup10pct,
                   depthdup20pct, depthdup50pct)

    dic["MAPPED_READS_DEDUP"] = dic_rmdup["[Total] Mapped Reads"]
    dic["TARGET_READS_DEDUP"] = dic_rmdup["[Target] Target Reads"]

    dic["DUPLICATE(%)"] = str2pctstr(
        1 - float(dic["MAPPED_READS_DEDUP"]) / float(dic["MAPPED_READS"]))
    dic["DUPLICATE_TARGET(%)"] = str2pctstr(
        1 - float(dic["TARGET_READS_DEDUP"]) / float(dic["TARGET_READS"]))

    dic["ON_TARGET_EXT(%)"] = dic_flank[
        "[Target] Fraction of Target Data in all data"]
    dic["ON_TARGET_READS_EXT(%)"] = dic_flank[
        "[Target] Fraction of Target Reads in all reads"]
    dic["CV_SCORE"] = calCV(region_sort, stat="mean")
    return dic


def bamdst_run(bamdstpath,
               sortbam,
               rmdupbam,
               bed,
               genome,
               bedtools,
               outdir="./",
               flank=100,
               mapQ=20,
               uncover=20):
    """Parallel run bamdst 3 times(sorted bam raw bed, sorted bam extend bed, rmdup bam raw bed)

    Arguments:
        bamdstpath {[str]} -- [path to bamdst]
        sortbam {[str]} -- [path to sorted bam]
        rmdupbam {[str]} -- [path to rndup bam]
        bed {[str]} -- [path to bed file]
        outdir {[str]} -- [path to temp dir]

    Keyword Arguments:
        flank {int} -- [base pair extend for bed file] (default: {100})
        mapQ {int} -- [description] (default: {20})
        uncover {int} -- [description] (default: {20})
    """

    rawdir = outdir + "/sort"
    flankdir = outdir + "/flank"
    rmdupdir = outdir + "/rmdup"
    flankbed = flankdir + "/flank.bed"
    makedir(rawdir, flankdir, rmdupdir)
    slopbed(bed, flank, flankbed, genome, bedtools)

    threads = []
    raw_process = threading.Thread(
        target=bamdst_subprocess,
        args=(
            bamdstpath,
            sortbam,
            bed,
            rawdir,
            mapQ,
            uncover,
        ))
    flank_process = threading.Thread(
        target=bamdst_subprocess,
        args=(
            bamdstpath,
            sortbam,
            flankbed,
            flankdir,
            mapQ,
            uncover,
        ))
    rmdup_process = threading.Thread(
        target=bamdst_subprocess,
        args=(
            bamdstpath,
            rmdupbam,
            bed,
            rmdupdir,
            mapQ,
            uncover,
        ))
    threads = [rmdup_process, raw_process, flank_process]
    startlog = [
        "Rmdup bam raw bed Start", "Sort bam raw bed bamdst Start",
        "Sort bam flank bed Start"
    ]
    finishlog = [
        "Rmdup bam raw bed Finished", "Sort bam raw bed Finished",
        "Sort bam flank bed Finished"
    ]
    for i in range(len(threads)):
        threads[i].setDaemon(True)
        threads[i].start()
        print ctime(), startlog[i]

    for i in range(len(threads)):
        threads[i].join()
        print ctime(), finishlog[i]

    sampleid = sortbam.split("/")[-1].split(".")[0]
    coverage_sort = outdir + "/sort/coverage.report"
    depth_distribution_sort = outdir + "/sort/depth_distribution.plot"
    insertsize_sort = outdir + "/sort/insertsize.plot"
    region_sort = outdir + "/sort/region.tsv.gz"
    coverage_rmdup = outdir + "/rmdup/coverage.report"
    depth_distribution_rmdup = outdir + "/rmdup/depth_distribution.plot"
    coverage_flank = outdir + "/flank/coverage.report"
    dic = bamdst_integrate(sampleid, coverage_sort, depth_distribution_sort,
                           insertsize_sort, region_sort, coverage_rmdup,
                           depth_distribution_rmdup, coverage_flank)

    return dic


def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s", "--sortedbam", help="sorted bam file", required=True)
    parser.add_argument(
        "-r", "--rmdupbam", help="rmdup bam file", required=True)
    parser.add_argument("-b", "--bed", help="target bed file", required=True)
    parser.add_argument(
        "-d",
        "--bamdst",
        help="path to bamdst software [/GPFS01/softwares/bamdst/bamdst]",
        default="/GPFS01/softwares/bamdst/bamdst")
    parser.add_argument(
        "-bd", "--bedtools", help="path to bedtools", required=True)
    parser.add_argument(
        "-g", "--genome", help="path to reference fai file", required=True)
    parser.add_argument(
        "-t", "--tempdir", help="path to bamdst output dir [./]", default="./")
    parser.add_argument(
        "-o", "--output", help="output file name [stdout]", default=sys.stdout)
    parser.add_argument(
        "-f",
        "--flank",
        help="base pairs in each direction [100]",
        default=100)
    parser.add_argument(
        "-m",
        "--mapq",
        help=
        "map quality cutoff value, greater or equal to the value will be count [20]",
        default=20)
    parser.add_argument(
        "-u",
        "--uncover",
        help="region will included in uncover file if below it [20]",
        default=20)
    parser.add_argument(
        "-n",
        "--noheader",
        action="store_true",
        help="print no header [False]")
    parser.add_argument(
        "-p", "--sep", help="delimiter, default:,", default=",")
    return parser.parse_args()


def headerlis():
    return [
        "#SAMPLE", "CLEAN_READS", "CLEAN_BASES", "INSERT_SIZE", "MAPQ20(%)",
        "DUPLICATE(%)", "DUPLICATE_TARGET(%)", "ON_TARGET_CORE(%)",
        "ON_TARGET_EXT(%)", "ON_TARGET_READS_CORE(%)",
        "ON_TARGET_READS_EXT(%)", "RATIO_OF_MAPPED(%)", "MEAN_DEPTH",
        "MEDIAN_DEPTH", "1X_COVERAGE(%)", "20X_COVERAGE(%)", "50X_COVERAGE(%)",
        "100X_COVERAGE(%)", "200X_COVERAGE(%)", "500X_COVERAGE(%)",
        "10%_COVERAGE(%)", "20%_COVERAGE(%)", "50%_COVERAGE(%)",
        "MEAN_DEPTH_DEDUP", "MEDIAN_DEPTH_DEDUP", "1X_COVERAGE_DEDUP(%)",
        "20X_COVERAGE_DEDUP(%)", "50X_COVERAGE_DEDUP(%)",
        "100X_COVERAGE_DEDUP(%)", "200X_COVERAGE_DEDUP(%)",
        "500X_COVERAGE_DEDUP(%)", "10%MEAN_COVERAGE_DEDUP(%)",
        "20%MEAN_COVERAGE_DEDUP(%)", "50%MEAN_COVERAGE_DEDUP(%)", "CV_SCORE"
    ]


def write_infos(dic, output=sys.stdout, noheader=False, sep=","):
    if sep == 't':
        sep = '\t'
    if output == sys.stdout:
        w = sys.stdout
    else:
        w = open(output, "w")

    if not noheader:
        w.write(sep.join(headerlis()) + "\n")
    lis = []
    for h in headerlis():
        lis.append(dic[h])
    w.write(sep.join(map(str, lis)) + "\n")


def main():
    if len(sys.argv) == 1:
        sys.argv.append("-h")
    args = parseargs()
    dic = bamdst_run(args.bamdst, args.sortedbam, args.rmdupbam, args.bed,
                     args.genome, args.bedtools, args.tempdir, args.flank,
                     args.mapq, args.uncover)
    write_infos(dic, args.output, args.noheader, args.sep)


if __name__ == '__main__':
    main()
