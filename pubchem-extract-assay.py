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


def process_bioassay_file(gzip_file, handles, annotations, sid2aid):
    
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
    drCols = {}
    acCol = None
    for i in assay['results']:
        tid = i['tid']
        if 'tc' in i and 'dr_id' in i['tc']:
            dr_id = i['tc']['dr_id']
            if tid not in drCols:
                drCols[tid] = {}
            drCols[tid] = (dr_id, i['tc']['concentration'])
        else:
            tcols[i['tid']] = i["name"]
        if 'ac' in i and i['ac']:
            acCol = tid
        
    for entry in data:
        adout = AssayData()
        adout.assay_id = aid
        if entry["sid"] in sid2aid:
            adout.compound_id = "CID:%s" % (sid2aid[entry["sid"]])
        else:
            adout.compound_id = "SID:%s" % (entry["sid"])
        drOut = {}
        if 'data' in entry:
            acVal = None
            if acCol is not None:
                for d in entry['data']:
                    tid = d["tid"]
                    if tid == acCol:
                        acVal = getValue(d["value"])
            for d in entry['data']:
                tid = d["tid"]
                if tid in tcols:
                    if 'fval' in d['value']:
                        adout.float_vals[tcols[tid]] = d['value']['fval']
                    if 'sval' in d['value']:
                        adout.string_vals[tcols[tid]] = d['value']['sval']
                    if 'ival' in d['value']:
                        adout.int_vals[tcols[tid]] = d['value']['ival']
                if tid in drCols:
                    dr_id, concentration = drCols[tid]
                    if dr_id not in drOut:
                        drOut[dr_id] = ResponseCurve()
                        drOut[dr_id].assay_id = aid
                        c = drOut[dr_id].compounds.add()
                        c.compound = adout.compound_id
                        c.ratio = 1.0
                        if acVal is not None:
                            s = drOut[dr_id].summary.add()
                            s.type = s.AC50
                            s.value = acVal
                    drV = drOut[dr_id].values.add()
                    drV.dose = concentration
                    drV.response = getValue(d["value"])
        if 'comment' in entry:
            adout.comment = entry['comment']
        #print aid, entry["sid"], entry["outcome"], e_data
        if "outcome" in entry and entry["outcome"] == "active":
            adout.active = True
        handles['AssayData'].write(json.dumps( json_format.MessageToDict(adout, including_default_value_fields=True)) + "\n")
        for dr in drOut.values():
            handles['AssayDoseResponse'].write(json.dumps( json_format.MessageToDict(dr, including_default_value_fields=True)) + "\n")
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dir", default=None)
    parser.add_argument("-f", "--gzip_file", default=None)
    parser.add_argument("-a", "--annotation", default=None, help="ftp://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/Extras/Aid2Annotation.gz")
    parser.add_argument("-s", "--sid2aid", default=None, help="ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SID.gz")
    args = parser.parse_args()
    
    annotations = {}
    if args.annotation is not None:
        with gzip.GzipFile(args.annotation) as handle:
            for line in handle:
                row = line.rstrip().split("\t")
                if row[0] not in annotations:
                    annotations[row[0]] = {}
                annotations[row[0]][row[1]] = row[2]
    sid2aid = {}
    if args.sid2aid is not None:
        with gzip.GzipFile(args.sid2aid) as handle:
            for line in handle:
                row = line.rstrip().split("\t")
                sid2aid[row[1]] = row[0]
        
    
    handles = {
        "Assay" : open("out.Assay.json", "w"),
        "AssayData" : open("out.AssayData.json", "w"),
        "AssayDoseResponse" : open("out.AssayDoseResponse.json", "w")
    }
    if args.gzip_file is not None:
        process_bioassay_file(args.gzip_file, handles, annotations, sid2aid)
    if args.dir is not None:
        for gzip_file in glob(os.path.join(args.dir, "*.json.gz")):
            print "Process %s" % (gzip_file)
            process_bioassay_file(gzip_file, handles, annotations, sid2aid)

    for i in handles.values():
        i.close()