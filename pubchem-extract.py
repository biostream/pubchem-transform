#!/usr/bin/env

import sys
from phenotype_pb2 import Compound
from urllib2 import urlopen
from google.protobuf import json_format

def message_to_json(message):
    return json.dumps(json_format.MessageToDict(message))

RECORD = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/record/JSON"
ASSAY_SUMMARY = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/assaysummary/JSON"
SYNONYMS = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/synonyms/JSON"

CID = sys.argv[1]
