#! /usr/bin/python
'''
Authors: Malisa Smith smimal@ohsu.edu, Ryan Spangler spanglry@ohsu.edu
Updated by: Theodore J. LaGrow lagrow@ohsu.edu

This program converts CTDD drug response information into
protobuf data based on the BMEG sample.proto schema.

Source: ftp://caftpd.nci.nih.gov/pub/OCG-DCC/CTD2/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/CTRPv2.0_2015_ctd2_ExpandedDataset.zip

The four files of interest (for this converter) are:
1) v20.data.curves_post_qc.txt
2) v20.meta.per_compound.txt
3) v20.meta.per_cell_line.txt
4) v20.meta.per_experiment.txt

Update by TL:
The data has been updated to incorporate the PubChem meta data. 
The additional file of interest (for this converter) is:
5) Drug_Conversion_MASTER_UPDATED.xlsx, sheet_name: Master_noerror

'''
########used interchangably
import phenotype_pb2
#from bmeg import phenotype_pb2 
#######

from google.protobuf import json_format
import json, sys, argparse, os
import csv #for drug data
import string
import re
import pandas
####### Added by TL
import xlrd #extracting data from the excel spreadsheet
import math #used checking for NaN
import time #testing efficiency 
#######

compound_gid = dict()

def parse_args(args):
    # We don't need the first argument, which is the program name
    args = args[1:]

    # Construct the parser
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # Now add all the options to it
    parser.add_argument('--response', type=str, help='Path to the drug response experiment data you want to import')
    parser.add_argument('--metadrug', type=str, help='Path to the drug meta data you want to import')
    parser.add_argument('--metacellline', type=str, help='Path to the cell line meta data you want to import')
    parser.add_argument('--metaexperiment', type=str, help='Path to the experiment meta data you want to import')
    parser.add_argument('--data', type=str, help='Path to the experiment data you want to import')
    ####### Added by TL
    parser.add_argument('--exceldata', type=str, help='Path to the Excel spreadsheet of data you want to import')
    parser.add_argument('--sheet_name', type=str, help='Name of sheet from Excel spreadsheet of data you want to import')
    #######
    parser.add_argument('--out', type=str, help='Path to output file (.json or .pbf_ext)')
    parser.add_argument('--multi', type=str, help='Path to output file (.json or .pbf_ext)')
    parser.add_argument('--format', type=str, default='json', help='Format of output: json or pbf (binary)')
    return parser.parse_args(args)

########################################

def find_biosample(state, source, barcode, sample_type):
    sample_name = 'biosample:CCLE:' + barcode
    biosample = state['Biosample'].get(sample_name)
    if biosample is None:
        biosample = schema.Biosample()
        biosample.name = sample_name
        biosample.dataset_id = "CCLE"
        biosample.source = source
        biosample.barcode = barcode
        biosample.sampleType = sample_type
        state['Biosample'][sample_name] = biosample

    return biosample

def append_unique(l, i):
    if not i in l:
        l.append(i)

def CheckUnicode(s):
    """ Checking if a string is unicode/NaN to not raise an error """
    try:
        math.isnan(s)
        return unicode('', "UTF-8")
    except:
        if isinstance(s, float):
            s = str(int(s))
        if isinstance(s, unicode):
            return s
        else:
            return unicode(s, "UTF-8")


def gid_naming(name, cid, sid):
    #print "Cid: ", cid
    if cid != unicode('', "UTF-8") and cid[:1] != "[":
        return "compound:CID" + cid, "compound"

    elif sid != unicode('', "UTF-8") and sid[:1] != "[":
        return "compound:SID" + sid, "substance"

    else:
        return "compound:" + name, "unknown"


def process_drugs(emit, input, excelData): #row is a namedtuple
    """ Modified by TL """
    compounds = set() #keeping track of unquie compounds
    m1 = dict() #for repetition of unique compounds when come across again
    m2 = dict()
    #sourceDict = {} #for testing

    start = time.time() #keeping track for testing

    for row in excelData.itertuples():
        if row.Name in m1:
            m1[row.Name] += 1
        else:
            m1[row.Name] = 1

    for key in m1:
        if m1[key] != 1:
            m2[key] = m1[key]

    between_row_time_s = time.time()
    times = set()
    counting_ctdd = 0

    for excelRow in excelData.itertuples():

        compound_name, pubchemtypedelineation = gid_naming(CheckUnicode(excelRow.Name), CheckUnicode(excelRow.PubChem_CID), CheckUnicode(excelRow.PubChem_SID))

        if compound_name not in compounds:
            compound_gid[excelRow.Name] = compound_name
            compound = phenotype_pb2.Compound()
            compound.gid = compound_name
            compound.id = compound_name
            compound.pubchemtype = pubchemtypedelineation
            compound.name = CheckUnicode(excelRow.Name)
            compound.smiles = CheckUnicode(excelRow.Smiles)
            compound.pubchemcid = CheckUnicode(excelRow.PubChem_CID)
            compound.pubchemsid = CheckUnicode(excelRow.PubChem_SID)
            compound.toxicity    = CheckUnicode(excelRow.Toxicity)
            compound.bioassays   = CheckUnicode(excelRow.BioAssays)
            compound.chebi_id    = CheckUnicode(excelRow.ChEBI_ID)
            compound.biological_role   = CheckUnicode(excelRow.Biological_Role)
            compound.biological_role_info   = CheckUnicode(excelRow.Biological_Role_Info)
            compound.application   = CheckUnicode(excelRow.Application)
            compound.application_info    = CheckUnicode(excelRow.Application_Info)
            compound.source = CheckUnicode(excelRow.Source) #changed back to single value, will be merged during ingestion
            compound.fda_approved = CheckUnicode(excelRow.fda_approved)
            compound.fda_approved_date = CheckUnicode(excelRow.fda_approved_date)
            compound.fda_data_origin = CheckUnicode(excelRow.fda_data_origin)
            compound.field_of_fda_approval_from = CheckUnicode(excelRow.field_of_fda_approval_from)

            syn = [excelRow.broad_cpd_id, excelRow.other_name, excelRow.GNF_REG_ID, excelRow.Drug_ID]
            for n in syn:
                if CheckUnicode(n) != unicode('', "UTF-8"): #chcek for no errors
                    compound.synonyms.append(CheckUnicode(n))

            #need to check one time and grab one time
            once = True
            current_compound_name = ""

            for row in input.itertuples():
                if once != True:
                    break
                if row.cpd_name == CheckUnicode(excelRow.Name) and once == True:
                    
                    current_compound_name = row.cpd_name
                    counting_ctdd += 1
                    e2 = time.time()
                    timestamp = e2-between_row_time_s
                    print "{}  {}, elapTime: {}, gid:{}".format(counting_ctdd, row.cpd_name, timestamp, compound_name)
                    times.add(timestamp)
                    between_row_time_s = time.time()


                    compound.status = row.cpd_status
                    target = row.gene_symbol_of_protein_target
                    if target and isinstance(target, str):
                        compound.target = str('gene:' + target)
                    compound.report = row.target_or_activity_of_compound
                    compound.rationale = row.inclusion_rationale

                    once = False
            input = input[input.cpd_name != current_compound_name] #this is to help speed up the process, should be log time complexity when parsing the dataframe


            emit(compound) #producing message

            if excelRow.Name in m2 and m2[excelRow.Name] != 1:
                m2[excelRow.Name] -= 1
            else:
                compounds.add(compound_name) #final step to add to the unique compound set

    print "Shortest: {}, Longest: {}, Average: {}".format(min(times), max(times), sum(times)/float(len(times)))


    #There shouldn't be anything left but just a final pass
    for row in input.itertuples():
        # create drug message for CTDD compound

        compound_name = "compound:" + CheckUnicode(row.cpd_name)
        

        if compound_name not in compounds:
            compound_gid[row.cpd_name] = compound_name
            compound = phenotype_pb2.Compound()
            compound.id = compound_name
            compound.gid = compound_name
            compound.name = row.cpd_name
            compound.smiles = row.cpd_smiles
            compound.status = row.cpd_status
            target = row.gene_symbol_of_protein_target
            if target and isinstance(target, str):
                compound.target = str('gene:' + target)
            compound.report = row.target_or_activity_of_compound
            compound.rationale = row.inclusion_rationale
            compound.synonyms.append(row.broad_cpd_id)

            emit(compound)
            compounds.add(compound_name)
            
    end = time.time()
    print "Processing time took: ", end-start

    

def process_response(emit, input, data):
    gid_set = set()
    ##### Added by TL, testing efficiency
    start = time.time()
    #####


    for row in input.itertuples():
        if isinstance(row.ccle_primary_site, str):
            site = row.ccle_primary_site.upper()
            sample_name = '%s_%s' % (row.ccl_name, site)
            sample = 'biosample:CCLE:' + sample_name
            ###### Added by TL
            compound_name = compound_gid[row.cpd_name]
            ######
            gid = "responseCurve:%s:%s" % (sample_name, row.cpd_name)

            if gid not in gid_set:
                response = phenotype_pb2.ResponseCurve()
                response.gid = gid
                response.responseType = phenotype_pb2.ResponseCurve.ACTIVITY
                response.compound = compound_name
                response.sample = sample

                s = response.summary.add()
                s.type = phenotype_pb2.ResponseSummary.EC50
                s.value = row.apparent_ec50_umol
                s.unit = "uM"
            
                s = response.summary.add()
                s.type = phenotype_pb2.ResponseSummary.AUC
                s.value = row.area_under_curve
                s.unit = "uM"

                for m in data.loc[lambda x: x.master_cpd_id==row.master_cpd_id, : ].loc[lambda x: x.experiment_id==row.experiment_id].itertuples():
                    dr = response.values.add()
                    dr.dose = m.cpd_conc_umol
                    dr.response = m.cpd_expt_avg_log2
            
                emit(response)
                gid_set.add(gid)
    ##### Added by TL
    end = time.time()
    print "Drug processing response curve took: ", end-start
    #####

def convert_all_ctdd(responsePath, metadrugPath, metacelllinePath, metaexperimentPath, dataPath, excelDataPath, sheetName, out, multi=None):
    """ Modified by TL """    
    # Read in Compound information into a pandas dataframe.
    compound_df = pandas.read_table(metadrugPath)
    # Read in Cell line information
    ccl_df = pandas.read_table(metacelllinePath)
    # Read in data curves for experiments
    datacurves_df = pandas.read_table(responsePath)
    # Read in meta experimental data
    metaexperiment_df = pandas.read_table(metaexperimentPath)
    
    
    print "\nReading excel data"
    data_xls = pandas.read_excel(excelDataPath, sheetName, index_col=None)
    print "  Converting to csv"
    data_xls.to_csv('intermcsvfile.csv', encoding='utf-8', index=False, sep="\t")
    print "    Finished converting excel to csv"


    exceldata_df = pandas.read_table('intermcsvfile.csv')

    print "      Merging tables"
    ctdd_merged = pandas.merge(datacurves_df, metaexperiment_df, how='left', on=['experiment_id']) # merge experiment data
    ctdd_merged = pandas.merge(ctdd_merged, compound_df, how='left', on=['master_cpd_id']) # merge with compound data frame
    #print "ctdd_merged: ", list(ctdd_merged.columns.values)
    #print "Excel:      ", list(exceldata_df.columns.values)
    

    #ctdd_merged = pandas.merge(ctdd_merged, exceldata_df, how='outer', on=['broad_cpd_id']) # merge all headers to overall df
    ctdd_merged = pandas.merge(ctdd_merged, ccl_df, how='left', on=['master_ccl_id']) # merge with cell line data frame
    #print "ctdd_merged final: ", list(ctdd_merged.columns.values)
    
    ctdd_data = pandas.read_table(dataPath)
    #print ctdd_merged
    
    out_handles = {}
    def emit_json_single(message):
        if 'main' not in out_handles:
            out_handles['main'] = open(out, "w")
        msg = json.loads(json_format.MessageToJson(message))
        msg["#label"] = message.DESCRIPTOR.full_name
        out_handles['main'].write(json.dumps(msg))
        out_handles['main'].write("\n")

    def emit_json_multi(message):
        if message.DESCRIPTOR.full_name not in out_handles:
            out_handles[message.DESCRIPTOR.full_name] = open(multi + "." + message.DESCRIPTOR.full_name + ".json", "w")
        msg = json.loads(json_format.MessageToJson(message))
        out_handles[message.DESCRIPTOR.full_name].write(json.dumps(msg))
        out_handles[message.DESCRIPTOR.full_name].write("\n")

    if out is not None:
        emit = emit_json_single
    if multi is not None:
        emit = emit_json_multi
    

    print "        Tables merged!"
    print "          Writing out"
    ctdd_merged.to_csv("test.out", sep="\t")   
    print "            Finished Writing out!"
 

    print "\nStarting processing drugs"
    process_drugs(emit, ctdd_merged, exceldata_df)
    print "  Finshed processing drugs"
    print "\n\nStarting processing response"
    #process_response(emit, ctdd_merged, ctdd_data)
    #print "  Finished processing response"
    
########################################

def message_to_json(message):
    msg = json.loads(json_format.MessageToJson(message))
    msg['#label'] = message.DESCRIPTOR.name
    return json.dumps(msg)


def convert_to_profobuf(responsePath, metadrugPath, metacelllinePath, metaexperimentPath, dataPath, excelDataPath, sheetName, out, multi):
    ####### Updated by TL, 7/7
    if responsePath and metadrugPath and metacelllinePath and metaexperimentPath and (out or multi) and format:
        convert_all_ctdd(responsePath, metadrugPath, metacelllinePath, metaexperimentPath, dataPath, excelDataPath, sheetName, out, multi)
    else:
        print("Please include all arguments")

    #write_messages(state, outpath, format)

if __name__ == '__main__':
    options = parse_args(sys.argv)
    ####### Updated by TL
    convert_to_profobuf(options.response, options.metadrug, options.metacellline, options.metaexperiment, dataPath=options.data, excelDataPath=options.exceldata, sheetName=options.sheet_name, out=options.out, multi=options.multi)
    ####### 

