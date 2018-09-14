#!/usr/bin/env python
import subprocess
import argparse
import sys


def gatk4_subproces(gatk4, bamfile, exacvcf, tempdir):
    """run gatk4 for calculate contamination
    
    Arguments:
        gatk4 {[str]} -- [path to gatk4]
        bamfile {[str]} -- [path to bamfle]
        exacvcf {[str]} -- [path to exacvcf]
        tempdir {[str]} -- [path to tempdir]
    
    Returns:
        [dict] -- [dict for results]
    """

    dic = {}
    tempfile = tempdir + "/getpileup.txt"
    output = tempdir + "/contamination.txt"
    cmdlis1 = [
        gatk4, "GetPileupSummaries", "-I", bamfile, "-V", exacvcf, "-O",
        tempfile
    ]
    cmdlis2 = [gatk4, "CalculateContamination", "-I", tempfile, "-O", output]
    subprocess.call(cmdlis1)
    subprocess.call(cmdlis2)
    sample = bamfile.split("/")[-1].split(".")[0]
    contamination = get_contamination(output)
    dic["#SAMPLE"] = sample
    dic["CONTAMINATION"] = contamination
    return dic


def get_contamination(file_contamination):
    """parse contamination output
    
    Arguments:
        file_contamination {[str]} -- [path to contamitation output]
    
    Returns:
        [float] -- [value of contamination]
    """

    fin = open(file_contamination).readlines()
    if len(fin) < 2:
        sys.stderr.write("Error: Contamination broken")
    contamination = round(float(fin[1].split()[1]), 4)
    return contamination


def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--bamfile", help="input bam file", required=True)
    parser.add_argument(
        "-v",
        "--exacvcf",
        help=
        "ExAC VCF file [/GPFS01/databases/GeneticDB/ExAC/release0.3.1/subsets/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz]",
        default=
        "/GPFS01/databases/GeneticDB/ExAC/release0.3.1/subsets/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz"
    )
    parser.add_argument(
        "-g",
        "--gatk4",
        help="path to gatk4 [/GPFS01/softwares/gatk-4.0.8.1/gatk]",
        default="/GPFS01/softwares/gatk-4.0.8.1/gatk")
    parser.add_argument(
        "-t", "--tempdir", help="path to bamdst output dir [./]", default="./")

    return parser.parse_args()


def main():
    if len(sys.argv) == 1:
        sys.argv.append("-h")
    args = parseargs()
    dic = gatk4_subproces(args.gatk4, args.bamfile, args.exacvcf, args.tempdir)
    print ",".join(map(str, [dic["#SAMPLE"], dic["CONTAMINATION"]]))


if __name__ == '__main__':
    main()