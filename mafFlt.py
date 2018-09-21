#!/usr/bin/env python
# File Name: mafFlt
# Author:Yong Wu
# Created Time: Fri 21 Sep 2018 02:22:15 PM CST
# Copyright: GENESEEQ Technology Inc. All rights reserved.
# Description:
#########################################################################
from collections import defaultdict
import os

def Ddict():
    return defaultdict(dict)

#文本文件转二维字典
def txt2dic(txt, headidx=0, keyidx=[1, 2, 3, 4]):
    dic = Ddict()
    header = open(txt).readlines()[headidx].rstrip().split("\t")
    for i in open(txt).readlines()[headidx + 1:]:
        t = i.rstrip().split("\t")
        key = tuple([t[x] for x in keyidx])
        for s in range(len(header)):
            dic[key][header[s]] = t[s]
    return dic, header

#maf转二维字典，可有version行或无version行
def maf2dic(inmaf):
    if open(inmaf).readline().startswith("#"):
        headidx = 1
    else:
        headidx = 0

    dic, header = txt2dic(inmaf, headidx, keyidx=[4, 5, 9, 11])

    for key, val in dic.items():
        t_vaf = float(val["t_alt_count"])/float(val["t_depth"])
        dic[key]["t_vaf"] = t_vaf
    return dic, header

#筛选函数
def filter_value(item_value, min, max):
    if min == None and max == None:
        return True
    if min == None:
        if item_value > max:
            return False
        else:
            return True
    elif max == None:
        if item_value < min:
            return False
        else:
            return True
    else:
        if min <= item_value <= max:
            return True
        else:
            return False


def str2float(s):
    if s in ["", "."]:
        return 0
    elif "%" in s:
        return float(s.replace("%", ""))
    else:
        return float(s)


#筛选maf某一项
def filter_item(inmaf, item, min=None, max=None):
    dic, header = maf2dic(inmaf)
    fail_dic = Ddict()
    for key in dic:
        value = str2float(dic[key][item])
        if filter_item(value, min, max):
            continue
        else:
            fail_dic[key] = dic[key]
    return fail_dic


#




"""
input file format[csv]:
maf_file_name
item1,min,max
item2,min,max
item3,min,max
item4,min,max
"""
