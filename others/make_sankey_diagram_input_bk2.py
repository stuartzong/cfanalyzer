#! /usr/bin/env python

import os, stat, os.path, time, datetime, subprocess
import re, sys, glob, argparse, csv, shutil, fileinput
from pprint import pprint
from itertools import islice
import itertools
import operator
import ConfigParser
from collections import defaultdict
import headers as HEADER


from jinja2 import Environment, FileSystemLoader
import logging
import colorlog

logger = colorlog.getLogger()
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
handler.setFormatter(
    colorlog.ColoredFormatter('%(log_color)s%(levelname)s:%(name)s:%(message)s'))
logger.addHandler(handler)


def get_files(bam_vcf_files, input_headers):
    """ 
    Dictionary: holding all file paths
    {patient ->
             {status ->
                     {file_identifier -> file_path}}}  
    """
    patient_files = dict()
    with open(bam_vcf_files, 'r') as fh:
        records = csv.DictReader(fh,  delimiter='\t')
        headers = set(records.fieldnames)
        input_headers = set(input_headers)
        # check if all mandatory columns prsent in input file
        if input_headers.issubset(headers):
            logger.info("input file has all mandatory columns! Continue...")
            for line in records:
                patient = line['patient']
                status = line['status']
                if patient not in patient_files:
                    patient_files[patient] = {}
                if status not in patient_files[patient]:
                    patient_files[patient][status] = line
                else:
                    logger.error('Duplicate status entries!'\
                                 'One tissue sequenced multiple times?')
                    sys.exit()
        else:
            logger.critical("Input file doesn't contain all mandatory headers!\n %s"
                            % list(input_headers))
            sys.exit()
        return patient_files

def make_sps_dict(sps_file, status):
    mutations = []
    sps_dict = dict()
    # primary_sps = "".join([wkdir, patient, "_primary.expands.sps"])
    # primary_snv= "".join([wkdir, patient, "_primary.expands.snv"])
    # refractory_sps = "".join([wkdir, patient, "_malignant.expands.sps" ])
    # refractory_snv = "".join([wkdir, patient, "_malignant.expands.snv" ])
    with open (sps_file, 'r') as fh:
        # skip the preheader line 
        next(fh)
        records = csv.DictReader(fh,  delimiter='\t')
        sps_list = []
        for line in records:
            chromosome = line["chr"]
            startpos = line["startpos"]
            subpopulation = line["SP"]
            mutation = "_".join([chromosome, startpos])
            # NA means EXPANDS failed to compute cellular frequency for this SNV
            if ("NA" not in subpopulation):
                sps_list.append(float(subpopulation))
                mutations.append(mutation)
                try:
                    sps_dict[mutation].append(subpopulation)
                except KeyError:
                    sps_dict[mutation] = [subpopulation]
        pprint(sps_dict)
        max_sps = max(sps_list)
        print "max_sps is:", max_sps
    return [sps_dict, max_sps, mutations]


def make_snv_dict(snv_file):
    snv_dict = dict()
    with open (snv_file, 'r' ) as fh:
     records = csv.DictReader(fh,  delimiter='\t')
     for line in records:
         chromosome = line["chr"]
         startpos = line["startpos"]
         gene = line["gene"]
         mutation = "_".join([chromosome, startpos])
         try:
             snv_dict[mutation].append(gene)
         except KeyError:
             snv_dict[mutation] = [gene]
    pprint(snv_dict)
    return snv_dict



def make_sankey_input_file(patient_files):
    for patient in patient_files:
        # sps = sub populations
        # primary_d = dict()
        # refractory_d = dict()
        # primary_snvd = dict()
        # refractory_snvd = dict()
        # wkdir = "/projects/trans_scratch/validations/workspace/szong/IF-AML/genome/expands/test/" 
        wkdir = os.getcwd()
        statuses = [i.lower() for i in patient_files[patient].keys()]
        print statuses
        if ('normal' in statuses and
            len(statuses) >= 3):
            # don't count normal 
            statuses = list(set(statuses)-set(['normal']))
            # interate status_combinations of statues
            # make sure to sort statuses in biopsy time point order 
            status_combinations = itertools.combinations(statuses, 2)
            for status_combination in status_combinations:
                print status_combination
                detail_dict = dict()
                for status in status_combination:
                    sps_dict = dict()
                    patient_status = '_'.join([patient, status])
                    directory = '/'.join([wkdir, patient_status])
                    sps_file = "".join([patient_status, '.expands.sps'])
                    sps_file = '/'.join([directory, sps_file])
                    snv_file = "".join([patient_status, '.expands.snv'])
                    # snv_file = '/'.join([directory, snv_file])
                    (sps_dict, max_sps, mutations) = make_sps_dict(sps_file, status)
                    snv_dict = make_snv_dict(snv_file)
                    detail_dict[status] = [sps_dict, max_sps, snv_dict, mutations]
                make_paired_sankey(patient, status_combination, detail_dict)

def make_paired_sankey(patient, status_combination, detail_dict):
    #adjust cellular frequency to 100% tumour content, write to sankey input file
    # sankey_header = ["primary_sp", "mutation", "refractory_sp", "value" ]
    patient_status = '_'.join(list(status_combination))
    sankey_file = ".".join([patient_status, "sankey" ])
    with open (sankey_file,  'wb') as fh:
        writer = csv.writer( fh, delimiter='\t' )
        status1 = status_combination[0]
        status2 = status_combination[1]
        sps1_dict = detail_dict[status1][0]
        sps2_dict = detail_dict[status2][0]
        max_sps1 = detail_dict[status1][1]
        max_sps2 = detail_dict[status2][1]
        snv1_dict = detail_dict[status1][2] 
        snv2_dict = detail_dict[status2][2] 
        mutations1 = detail_dict[status1][3]
        mutations2 = detail_dict[status2][3]
        mutations = list(set(mutations1 + mutations2))
        sps1_header = status1 + "_sp"
        sps2_header = status2 + "_sp"
        sankey_header = [sps1_header, 'mutation', sps2_header, 'value']
        writer.writerow(sankey_header)
        for mutation in mutations:
            # adjust cell frequency for status1
            try:
                sp1 = sps1_dict[mutation][0]
                if ("NA" not in sp1):
                    sp1 = "{:.3f}".format(float(sp1)/max_sps1)
                    sp1 = "_".join([status1,str(sp1)])
            except KeyError:
                sp1 = '_'.join([status1, 'NA'])
            # adjust cell frequency for status2
            try:
                sp2 = sps2_dict[mutation][0]
                if ("NA" not in sp2):
                    sp2 = "{:.3f}".format(float(sp2)/max_sps2)
                    sp2 = "_".join([status2,str(sp2)])
            except KeyError:
                sp2 = '_'.join([status2, 'NA'])
            # get gene names from snv dictionary    
            try:
                gene = snv1_dict[mutation][0]
            except KeyError:
                gene = snv2_dict[mutation][0]
            mutation = "_".join([gene, mutation])
            sankey_width = 0.075
            if ("MODIFIER" not in mutation):
                sankey_width = 0.2 
            writer.writerow([sp1, mutation, sp2, sankey_width] )


def modify_sankey(old_sankey, new_sankey):
    d = dict()
    with open (new_sankey,  'wb') as fh:
        writer = csv.writer( fh, delimiter='\t' )
        with open (old_sankey, 'r') as handle:
             records = csv.DictReader(handle,  delimiter='\t')
             headers = records.fieldnames
             writer.writerow(headers)
             for line in records:
                 primary_sp = line['primary_sp']
                 mutation = line["mutation"]
                 gene = mutation.split('_')[0]
                 refractory_sp = line["refractory_sp"]
                 value = line["value"]
                 if ('MODIFIER' not in mutation):
                     #content = [line[i] for i in headers]
                     content = [primary_sp, gene, refractory_sp, value]
                     print "This is a coding mutation!"
                     writer.writerow(content)
                 else: #make a modifier dictionary
                     key = ':'.join([primary_sp, refractory_sp])
                     try:
                         d[key].append(mutation)
                     except KeyError:
                         d[key] = [mutation]
        pprint(d)
        #write the number of modifers in each group into the sankey file
        for key in d:
            sl_key = key.split(':')
            primary_sp = sl_key[0]
            primary_cf = primary_sp.split('_')[1]
            refractory_sp = sl_key[1]
            refractory_cf = refractory_sp.split('_')[1]
            num_modifiers = len(list(set(d[key])))
            mutation = '_'.join([str(num_modifiers), 'MODIFIERS', primary_cf, refractory_cf])
            value = int(num_modifiers) * 0.02
            content = [primary_sp, mutation, refractory_sp, str(value)]
            writer.writerow(content)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Make sankey diagram input files!')
    parser.add_argument(
        '-i', '--meta_file',
        help='specify meta file, which gives patient, status info',
        required=True)
    args = parser.parse_args()
    return args


def __main__():
    print "Generating patient_files dictionary!"
    args = parse_args()
    meta_file = args.meta_file
    #input_files = "bam_vcf_test.tmp"
    input_headers = HEADER.INPUT_FILE_HEADER
    patient_files = get_files(meta_file, input_headers)
    pprint(patient_files)
    # patient_files = make_patient_vcf_file_dict(meta_file)

    #variant_summary = "SNV_summary_with_normal.txt.filtered"
    #variant_summary = "SNV_summary_with_normal.txt.all"
    make_sankey_input_file(patient_files)

if __name__ == '__main__':
    __main__()


