#!/usr/bin/env python

import argparse
import json
import gzip
from phenotype_pb2 import Compound
from urllib2 import urlopen
from urllib import quote
from google.protobuf import json_format
#from zeep import Client

def message_to_json(message):
    return json.dumps(json_format.MessageToDict(message))

RECORD = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/record/JSON"
ASSAY_SUMMARY = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/assaysummary/JSON"
SYNONYMS = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/synonyms/JSON"
ASSAY_RECORD = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/%s/record/JSON"
DESCRIPTION = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/description/JSON"
PUBMED = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/xrefs/pubmedid/JSON"
NAME_CID = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/cids/JSON"

def run_search(args):
    u = NAME_CID % (quote(args.name))
    try:
        record_txt = urlopen(u).read()
    except:
        return
    record = json.loads(record_txt)
    record['Name'] = args.name
    print json.dumps(record)

def run_extract(args):
    """
    #ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz
    GENE_FILE = sys.argv[1]
    gene_map = {}
    with gzip.GzipFile(GENE_FILE) as handle:
        for line in handle:
            row = line.rstrip().split("\t")
            gene_map[row[1]] = row[2]

    CID = sys.argv[2]

    record_txt = urlopen(RECORD % CID).read()
    print record_txt
    """
    for CID in args.ids:
        record_txt = urlopen(RECORD % CID).read()
        record = json.loads(record_txt)
        for compound in record["PC_Compounds"]:
            #print compound.keys()
            props = {}
            for prop in compound["props"]:
                if "name" in prop["urn"]:
                    label = prop["urn"]["label"] + ":" +  prop["urn"]["name"]
                else:
                    label = prop["urn"]["label"]
                for i in prop["value"].values():
                    props[label] = i
            print props

        syn_txt = urlopen(SYNONYMS % CID).read()
        syn = json.loads(syn_txt)
        print syn

        assay_txt = None
        try:
            assay_txt = urlopen(ASSAY_SUMMARY % CID).read()
        except:
            pass

        if assay_txt is not None:
            assay = json.loads(assay_txt)

            columns = assay["Table"]["Columns"]["Column"]

            for row in assay["Table"]["Row"]:
                cells = row["Cell"]
                data = dict(zip(columns, cells))
                if len(data['Target GeneID']):
                    print data

    #chebi_id = sys.argv[1]
    #client = Client('http://www.ebi.ac.uk/webservices/chebi/2.0/webservice?wsdl')
    #print client.service.getCompleteEntity(chebi_id)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    subparser = parser.add_subparsers()

    parser_search = subparser.add_parser("search")
    parser_search.add_argument("name")
    parser_search.set_defaults(func=run_search)

    parser_sync = subparser.add_parser("extract")
    parser_sync.add_argument("ids", nargs="+")
    parser_sync.set_defaults(func=run_extract)


    args = parser.parse_args()
    args.func(args)
