#!/usr/bin/python -u
# -*-coding:utf-8 -*-
__author__ = 'wuy'

import re
import sys
import time

import requests

import xmltodict

"""[animalModel
molecularGenetics
heterogeneity]

Returns:
    [type] -- [description]
"""

s = requests.session()
s.keep_alive = False
requests.adapters.DEFAULT_RETRIES = 100

reload(sys)
sys.setdefaultencoding('utf-8')
headers = {
    "User-Agent":
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:57.0) Gecko/20100101 Firefox/57.0",
    "Accept":
    "application/json, text/javascript, */*; q=0.01",
    "Accept-Encoding":
    "gzip, deflate, br",
    "ApiKey":
    "G13AC0SAQ_Se7rYUERJVpg"
}
cookies = {
    "ApiKey": "G13AC0SAQ_Se7rYUERJVpg",
    "_ga": "GA1.2.316904006.1516673540",
    "_gid": "GA1.2.879590851.1517193014"
}


def clinicalSynopsis(mim):
    url = "http://api.omim.org/api/clinicalSynopsis?mimNumber=%s&include=clinicalSynopsis&format=xml" % mim
    try:
        r = requests.get(url, headers=headers, cookies=cookies)
        #return r.json()["omim"]["clinicalSynopsisList"][0]["clinicalSynopsis"]
        dic = xmltodict.parse(
            r.text)["omim"]["clinicalSynopsisList"]["clinicalSynopsis"]
        #print dic
        lis = []
        for key, val in dic.items():
            t = val.split("\n")
            for i in t:
                lis.append(",".join([key, re.sub("{.+}", "", i)]))
        return "###".join(lis)
    except:
        return "NA"


def test(mim):
    url = "http://api.omim.org/api/clinicalSynopsis?mimNumber=%s&include=clinicalSynopsis&format=xml" % mim
    try:
        if requests.get(url, headers=headers, cookies=cookies):
            return True
    except:
        return False


def wait(s):
    if test("100100"):
        time.sleep(10)
    else:
        sys.stderr.write("sleep\n")
        time.sleep(s)
        wait(s)


# first time run
def main(mim2gene):
    with open(mim2gene) as fin:
        for i in fin:
            if i.startswith("#"):
                continue
            t = i.split("\t")
            if "phenotype" in t[1]:
                wait(3600)
                clinical = clinicalSynopsis(t[0])
                print "\t".join([t[0], clinical])
            else:
                continue



"""
if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print "py mim2gene.txt"
        exit(1)
    main(sys.argv[1])
"""
