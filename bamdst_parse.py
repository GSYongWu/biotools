#!/usr/bin/env python
import os
import subprocess
import threading

from pybedtools import BedTool


def makedir(*dirs):
    for path in dirs:
        if not os.path.exists(path):
            os.makedirs(path)


def str2pctstr(floatstr):
    return str(round(float(floatstr) * 100, 2)) + "%"


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


def bamdst_subprocess(bamdrspath, bam, bed, outdir, mapQ=20, uncover=20):
    """Run bamdrs command in subprocess 
    
    Arguments:
        bamdrspath {[string]} -- [path to bamdrs]
        bam {[string]} -- [path to bam file]
        bed {[string]} -- [path to bed file]
        outdir {[string]} -- [path to outdir]
    
    Keyword Arguments:
        mapQ {int} -- [map qual] (default: {20})
        uncover {int} -- [description] (default: {20})
    """

    cmdlis = [
        bamdrspath, '-q',
        str(mapQ), '--uncover',
        str(uncover), '-p', bed, '-o', outdir, bam
    ]
    subprocess.call(cmdlis)


def bamdst_run(bamdrspath,
               sortbam,
               rmdupbam,
               bed,
               outdir,
               flank=100,
               mapQ=20,
               uncover=20):
    """Parallel run bamdrs 3 times(sorted bam raw bed, sorted bam extend bed, rmdup bam raw bed)
    
    Arguments:
        bamdrspath {[string]} -- [path to bamdrs]
        sortbam {[string]} -- [path to sorted bam]
        rmdupbam {[string]} -- [path to rndup bam]
        bed {[string]} -- [path to bed file]
        outdir {[string]} -- [path to temp dir]
    
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
    slopbed(bed, flank, flankbed)

    threads = []
    raw_process = threading.Thread(
        target=bamdst_subprocess,
        args=(
            bamdrspath,
            sortbam,
            bed,
            rawdir,
            mapQ,
            uncover,
        ))
    flank_process = threading.Thread(
        target=bamdst_subprocess,
        args=(
            bamdrspath,
            sortbam,
            flankbed,
            flankdir,
            mapQ,
            uncover,
        ))
    rmdup_process = threading.Thread(
        target=bamdst_subprocess,
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
        #t.setDaemon(True)
        t.start()
    for t in threads:
        t.join()


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


def nxcoverage(depth_distribution, *ndepth):
    """Calculate coverage when depth >=N and median depth.
    
    Arguments:
        depth_distribution {string} -- [bamdst output depth_distribution.plot file]
    
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


def bamdst_integrate(sampleid, coverage_sort, depth_distribution_sort,
                     insertsize_sort, coverage_rmdup, depth_distribution_rmdup,
                     coverage_flank):

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
        float(dic["MAPPED_READS_DEDUP"]) / float(dic["MAPPED_READS"]))
    dic["DUPLICATE_TARGET(%)"] = str2pctstr(
        float(dic["TARGET_READS_DEDUP"]) / float(dic["TARGET_READS"]))

    dic["ON_TARGET_EXT(%)"] = dic_flank[
        "[Target] Fraction of Target Data in all data"]
    dic["ON_TARGET_READS_EXT(%)"] = dic_flank[
        "[Target] Fraction of Target Reads in all reads"]
    return dic


#[bamdrs_run and bamdrs_integrate test]


def header():
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
        "20%MEAN_COVERAGE_DEDUP(%)", "50%MEAN_COVERAGE_DEDUP(%)"
    ]


def main():
    sortbam = "/GPFS04/GSPipeline4/20180901_E00516_0382_AHLK35CCXY/HQData/Sample_B180621127568-KY400-2/B180621127568-KY400-2.sorted.bam"
    rmdupbam = "/GPFS04/GSPipeline4/20180901_E00516_0382_AHLK35CCXY/HQData/Sample_B180621127568-KY400-2/B180621127568-KY400-2.sorted.rmdup.bam"
    bamdst = "/GPFS01/softwares/bamdst/bamdst"
    bedfile = "/GPFS01/databases/GSCAP/db_for_201801_425/Selected_201801-425.raw.bed"

    sampleid = "B180621127568-KY400-2"
    coverage_sort = "sort/coverage.report"
    depth_distribution_sort = "sort/depth_distribution.plot"
    insertsize_sort = "sort/insertsize.plot"
    coverage_rmdup = "rmdup/coverage.report"
    depth_distribution_rmdup = "rmdup/depth_distribution.plot"
    coverage_flank = "flank/coverage.report"

    bamdst_run(bamdst, sortbam, rmdupbam, bedfile, "./")
    dic = bamdst_integrate(sampleid, coverage_sort, depth_distribution_sort,
                           insertsize_sort, coverage_rmdup,
                           depth_distribution_rmdup, coverage_flank)
    head = header()
    print ",".join(head)
    lis = []
    for h in head:
        lis.append(dic[h])
    print ','.join(map(str, lis))


if __name__ == '__main__':
    main()
