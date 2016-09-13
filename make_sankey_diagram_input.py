#! /usr/bin/env python

import os, stat, os.path, time, datetime, subprocess
import re, sys, glob, argparse, csv, shutil, fileinput
from pprint import pprint
# from itertools import islice
import itertools
# import operator
# import ConfigParser
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
            logger.critical("Input file doesn't contain all mandatory headers! %s"
                            % list(input_headers))
            sys.exit()
        return patient_files

def make_sps_dict(sps_file, status):
    mutations = []
    sps_dict = dict()
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
        max_sps = max(sps_list)
    return [sps_dict, max_sps, mutations]


def make_snv_dict(snv_file):
    snv_dict = dict()
    with open(snv_file, 'r') as fh:
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
    return snv_dict


def make_sankey_input_file(patient_files, template_dir, biopsy_time_points):
    sankey_files = []
    for patient in patient_files:
        wkdir = os.getcwd()
        statuses = [i.lower() for i in patient_files[patient].keys()]
        # print statuses
        if ('normal' in statuses and len(statuses) >= 3):
            # don't count normal 
            statuses = list(set(statuses)-set(['normal']))
            # iterate status_combinations of statues
            # make sure to sort statuses in biopsy time point order 
            # manually sort time point here 
            # statuses = ['diagnosis', 'biop1', 'biop2']
            statuses = biopsy_time_points
            status_combinations = itertools.combinations(statuses, 2)
            for status_combination in status_combinations:
                # print status_combination
                detail_dict = dict()
                for status in status_combination:
                    sps_dict = dict()
                    patient_status = '_'.join([patient, status])
                    directory = '/'.join([wkdir, patient_status])
                    sps_file = "".join([patient_status, '.expands.sps'])
                    sps_file = '/'.join([directory, sps_file])
                    snv_file = "".join([patient_status, '.expands.snv'])
                    (sps_dict, max_sps,
                     mutations) = make_sps_dict(sps_file, status)
                    snv_dict = make_snv_dict(snv_file)
                    detail_dict[status] = [sps_dict, max_sps,
                                           snv_dict, mutations]
                sankey_file = make_paired_sankey(patient, status_combination,
                                                 detail_dict, wkdir, template_dir)
                sankey_files.append(sankey_file)
    return sankey_files


def make_paired_sankey(patient, status_combination,
                       detail_dict, wkdir, template_dir):
    # adjust cellular frequency to 100% tumour content 
    patient_status = '_'.join([patient] + list(status_combination))
    sankey_file = ".".join([patient_status, "sankey"])
    modified_sankey = ".".join([patient_status, "sankey.modified"])
    with open(sankey_file,  'wb') as fh:
        writer = csv.writer(fh, delimiter='\t')
        (status1, status2) = status_combination
        sps1_header = status1 + "_sp"
        sps2_header = status2 + "_sp"
        (sps1_dict, max_sps1, snv1_dict, mutations1) = detail_dict[status1]
        (sps2_dict, max_sps2, snv2_dict, mutations2) = detail_dict[status2]
        mutations = list(set(mutations1 + mutations2))
        sankey_header = [sps1_header, 'mutation', sps2_header, 'value']
        writer.writerow(sankey_header)
        # make sankey_diagram.r script
        populate_template(patient_status, wkdir, modified_sankey,
                          sps1_header, template_dir)
        for mutation in mutations:
            subpopulation1 = adjust_cell_frequency(sps1_dict,
                                                   mutation,
                                                   max_sps1,
                                                   status1)
            subpopulation2 = adjust_cell_frequency(sps2_dict,
                                                   mutation,
                                                   max_sps2,
                                                   status2)
            gene = get_gene_name(snv1_dict, snv2_dict, mutation)
            mutation = "_".join([gene, mutation])
            sankey_width = 0.075
            if ("MODIFIER" not in mutation):
                sankey_width = 0.2 
            writer.writerow([subpopulation1,
                             mutation,
                             subpopulation2,
                             sankey_width])
    return sankey_file


def get_gene_name(snv1_dict, snv2_dict, mutation):
    try:
        gene = snv1_dict[mutation][0]
    except KeyError:
        gene = snv2_dict[mutation][0]
    return gene


def adjust_cell_frequency(subpopulation_dict, mutation, max_sps, status):
    try:
        subpopulation = subpopulation_dict[mutation][0]
        if ("NA" not in subpopulation):
            subpopulation = "{:.2f}".format(float(subpopulation)/max_sps)
            subpopulation = "_".join([status, str(subpopulation)])
    except KeyError:
        subpopulation = '_'.join([status, 'NA'])
    return subpopulation


def populate_template(patient_status, wkdir, sankey_input,
                      status1sp, template_dir):
    sankey_script = ".".join([patient_status, "sankey.r"])
    jinja2_env = Environment(loader=FileSystemLoader([template_dir]),
                             trim_blocks=True)
    template = jinja2_env.get_template('sankey_diagram_template.r')
    with open(sankey_script, 'wb') as opf:
        content = template.render(wkdir=wkdir,
                                  sankey_input=sankey_input,
                                  status1sp=status1sp)
        opf.write(content)
        logger.info('templated {0}'.format(sankey_script))
    return sankey_script


def modify_sankey(old_sankey, new_sankey):
    # count number of mutations for each group for better visualization 
    cf_combination_dict = dict()
    with open(new_sankey,  'wb') as fh:
        writer = csv.writer(fh, delimiter='\t')
        with open(old_sankey, 'r') as handle:
            records = csv.DictReader(handle,  delimiter='\t')
            headers = records.fieldnames
            writer.writerow(headers)
            for line in records:
                sps1 = line[headers[0]]
                sps2 = line[headers[2]]
                mutation = line[headers[1]]
                # gene = mutation.split('_')[0]
                value = line[headers[3]]
                if ('MODIFIER' in mutation):
                    pass
                else:
                    # dict to count mutations for each CF combinations
                    key = ':'.join([sps1, sps2])
                    try:
                        cf_combination_dict[key].append(mutation)
                    except KeyError:
                        cf_combination_dict[key] = [mutation]
        # write the number of mutatons in each group into the sankey file
        for key in cf_combination_dict:
            sl_key = key.split(':')
            sps1 = sl_key[0]
            sps2 = sl_key[1]
            sps1_cf = sps1.split('_')[1]
            sps2_cf = sps2.split('_')[1]
            num_mutations = len(list(set(cf_combination_dict[key])))
            mutation = '_'.join([str(num_mutations),
                                 'HM_muatations',
                                 sps1_cf,
                                 sps2_cf])
            value = int(num_mutations) * 0.05
            content = [sps1, mutation, sps2, str(value)]
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
    template_dir = '/home/szong/projects/development/expands/'
    input_headers = HEADER.INPUT_FILE_HEADER
    patient_files = get_files(meta_file, input_headers)
    pprint(patient_files)
    biopsy_time_points = ['diagnosis', 'biop1', 'biop2']
    sankey_files = make_sankey_input_file(patient_files,
                                          template_dir,
                                          biopsy_time_points)
    for sankey_file in sankey_files:
        new_sankey = '.'.join([sankey_file, 'modified'])
        modify_sankey(sankey_file,  new_sankey)


if __name__ == '__main__':
    __main__()
 

