import argparse
import sys

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--fastqinput", help="fastq file path list file", required=True)
    parser.add_argument("-sc","--config", help="config file", required=True)
    #parser.add_argument("-k", "--kit", help="kit infos file", required=True)
    #parser.add_argument("-c", "--crosslink", help="crosslink file", required=True)
    #parser.add_argument("-tl", "--tilecheck", help="tile check file", required=True)
    parser.add_argument(
        "-s", "--sortedbam", help="sorted bam file", required=True)
    parser.add_argument(
        "-r", "--rmdupbam", help="rmdup bam file", required=True)
    parser.add_argument("-b", "--bed", help="target bed file", required=True)
    
    parser.add_argument(
        "-t", "--tempdir", help="path to bamdst output dir [./]", default="./")
    parser.add_argument(
        "-o", "--output", help="output file name [stdout]", default=sys.stdout)
    
    return parser.parse_args()


