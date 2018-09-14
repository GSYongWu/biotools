#!/usr/bin/env python
import json
import os
import re
import sys

import numpy as np


def str2float(s):
    if s in ["", ".", "-"]:
        return None
    elif "%" in s:
        return float(s.replace("%", ""))
    else:
        return float(s)


# sample type
def sampletype(sampleid):
    dic = {}
    dic["TYPE"] = ""
    sampleid = ""
    if sampleid.startswith("P") or sampleid.startswith(
            "HP") or sampleid.startswith("ct") or sampleid.startswith(
                "MP") or (sampleid.startswith("D")
                          and sampleid.endswith("KY092-T")):
        dic["TYPE"] = "ctDNA"
    elif sampleid.startswith("ST") or sampleid.startswith("BT"):
        dic["TYPE"] = "Tumor"
    elif sampleid.startswith("HB") or sampleid.endswith(
            "MB") or sampleid.endswith("YB") or sampleid.endswith(
                "NA") or sampleid.startswith("B"):
        dic["TYPE"] = "Normal"
    else:
        dic["TYPE"] = "Tumor"
    return dic


# fastq site and lanes
def raw_fastq(input_fastqs):
    list_fastq = []
    list_lane = []
    with open(input_fastqs) as fin:
        for i in fin:
            i = i.rstrip()
            list_fastq.append(i)
            filename = i.split("/")[-1]

            if re.search("_(L00\d)_", filename):
                lane = re.search("_(L00\d)_", filename).group(1)
                if lane in list_lane:
                    list_lane.append(lane)
            else:
                lane = "UNK"
                if lane in list_lane:
                    list_lane.append(lane)
    return list_fastq, list_lane


# FC infos
def FCinfos(lane, kit, root="/GPFS01/SequencingData", machine="hiseq"):
    josnfile = os.path.join(root, kit, "BCLTMP", lane,
                            "Data/Intensities/BaseCalls/Stats/Stats.json")
    if os.path.exists(josnfile):
        jsondata = json.loads(open(josnfile)).read()
        return jsondata["ConversionResults"][0]["TotalClustersPF"]
    else:
        sys.stderr.write("Error: json file not exist!!")


def csv2nparray(csvfile):
    lis = []
    with open(csvfile) as csv:
        for i in csv:
            i = i.rstrip()
            if i.startswith("#"):
                continue
            t = i.split(",")
            lis.append(t)
    return np.array(lis)


# Raw fastq infos
def combineFastqInfos(infoarray):
    """Combine Raw Fastq Infos
    
    Arguments:
        infoarray {[np.array]} -- [np.array of raw infos]
    
    Returns:
        [dict] -- [infos for raw fastq]
    """

    dic = {}
    dic["PF_READS"] = round(np.sum(infoarray[:, 1].astype(np.float64)) / 2, 2)
    dic["LENGTH"] =  round(np.average(infoarray[:, 2].astype(np.float64)), 2)
    dic["PF_BASES"] =  round(np.sum(infoarray[:, 3].astype(np.float64)), 2)
    dic["GC(%)"] =  round(np.average(infoarray[:, 4].astype(np.float64)), 2)
    dic["Q20(%)"] =  round(np.average(infoarray[:, 5].astype(np.float64)), 2)
    dic["Q30(%)"] =  round(np.average(infoarray[:, 6].astype(np.float64)), 2)
    return dic


# crosslink
def crosslink(chean_bam_csv):
    dic = {}
    dic["CROSSLINK(%)"] = open(chean_bam_csv).readlines()[-1].rstrip().split(
        ",")[-1] + "%"
    return dic


# lims id
def get_limsid(sampleid, samplesheet):
    dic_ss = {}
    dic = {}
    with open(samplesheet) as ss:
        for i in ss:
            i = i.rstrip()
            if re.match("\d,", i):
                t = i.split(",")
                dic_ss[t[2]] = t[-1]
    if sampleid in dic_ss:
        dic["LIMS_ID"] = dic_ss[sampleid]
    else:
        dic["LIMS_ID"] = "."
    return dic


# stat
def stat(value, warning, item, threadtype="min", fail=None):
    """stat QC item PASS or Warning or Failed
    
    Arguments:
        value {[float]} -- [description]
        warning {[float or list]} -- [description]
        item {[str]} -- [description]
    
    Keyword Arguments:
        threadtype {str} -- [[min, mid, max]] (default: {"min"})
        fail {[float or list]} -- [description] (default: {None})
    
    Returns:
        [str] -- [return PASS or Warning or Failed]
    """

    if fail == "" or fail == 0 or fail == "NA":
        fail = None
    # if value is None, return PASS
    if value in ["", "-", "NA", None]:
        return "PASS"
    if warning in ["", "-", "NA", None, 0]:
        return "PASS"

    if threadtype == "min":
        if fail and value < fail:
            return "%s Failed" % item
        elif value < warning:
            return "%s Warning" % item
        else:
            return "PASS"
    elif threadtype == "mid":  # fail and warning is a list
        if fail and (value > fail[1] or value < fail[0]):
            return "%s Failed" % item
        elif value > warning[1] or value < warning[0]:
            return "%s Warning" % item
        else:
            return "PASS"
    else:
        if fail and value > fail:
            return "%s Failed" % item
        elif value > warning:
            return "%s Warning" % item
        else:
            return "PASS"
