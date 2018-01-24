#!/usr/bin/env python

import os
import re
import argparse
import json
import gzip
from glob import glob
from phenotype_pb2 import Compound, Assay, ResponseCurve, AssayData

from google.protobuf import json_format



def getValue(i):
    for d in ["sval", "fval", "ival"]:
        if d in i:
            return i[d]
    return None


def process_bioassay_file(gzip_file, handles, annotations):
    
    with gzip.GzipFile(gzip_file) as handle:
        data = json.loads(handle.read())
    
    entry = data['PC_AssaySubmit']
    
    assay = entry['assay']['descr']
    data = entry['data']
    #print json.dumps(assay, indent=4)
    #print json.dumps(data, indent=4)
    
    aid = str(assay["aid"]["id"])

    aout = Assay()
    aout.id = aid
    aout.name = assay["name"]
    if "xref" in assay:
        for x in assay["xref"]:
            if "pmid" in x["xref"]:
                aout.pubmed.append(str(x["xref"]["pmid"]))
            if "gene" in x["xref"]:
                aout.gene_target.append(str(x["xref"]["gene"]))                
    if aid in annotations:
        if "Assay Cell Type" in annotations[aid]:
            aout.cellline = annotations[aid]["Assay Cell Type"]
    handles['Assay'].write(json.dumps( json_format.MessageToDict(aout)) + "\n")
    handles['Assay'].flush()
    
    tcols = {}
    for i in assay['results']:
        tcols[i['tid']] = i["name"]

    e_data = {}
    for entry in data:
        adout = AssayData()
        if 'data' in entry:
            for d in entry['data']:
                e_data[ tcols[d["tid"]] ] = getValue(d["value"])
        if 'comment' in entry:
            adout.comment = entry['comment']
        #print aid, entry["sid"], entry["outcome"], e_data
        adout.assay_id = aid
        adout.compound_id = "SID:%s" % (entry["sid"])
        if "outcome" in entry and entry["outcome"] == "active":
            adout.active = True
        handles['AssayData'].write(json.dumps( json_format.MessageToDict(adout, including_default_value_fields=True)) + "\n")

    
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dir", default=None)
    parser.add_argument("-f", "--gzip_file", default=None)
    parser.add_argument("-a", "--annotation", default=None)
    
    args = parser.parse_args()
    
    annotations = {}
    if args.annotation is not None:
        with gzip.GzipFile(args.annotation) as handle:
            for line in handle:
                row = line.rstrip().split("\t")
                if row[0] not in annotations:
                    annotations[row[0]] = {}
                annotations[row[0]][row[1]] = row[2]
    
    handles = {
        "Assay" : open("out.Assay.json", "w"),
        "AssayData" : open("out.AssayData.json", "w"),
    }
    if args.gzip_file is not None:
        process_bioassay_file(args.gzip_file, handles, annotations)
    if args.dir is not None:
        for gzip_file in glob(os.path.join(args.dir, "*.json.gz")):
            print "Process %s" % (gzip_file)
            process_bioassay_file(gzip_file, handles, annotations)

    for i in handles.values():
        i.close()