#! /usr/bin/env python

"""
This script does the following:
-> make expands input files
-> run expands
"""


import os, stat, os.path, time, datetime, subprocess
import re, sys, glob, argparse, csv, shutil, fileinput
from pprint import pprint
from itertools import islice
import operator
import ConfigParser
from collections import defaultdict


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


def get_files(bam_vcf_files):
    """ 
    Dictionary: holding all file paths
    {patient ->
             {status ->
                     {file_identifier -> file_path}}}  
    """
    patient_files = dict()
    with open(bam_vcf_files, 'r') as fh:
        records = csv.DictReader(fh,  delimiter='\t')
        # headers = records.fieldname
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
        return patient_files


def example_display(n, iterable):
    """ 
    Return first n items of the iterable as a list! 
    """
    return list(islice(iterable, n))


def group_readable(file_path):
    m = os.stat(file_path).st_mode
    #other_execute  = bool(m & 0001)
    #other_write = bool(m & 0002)
    #other_read  = bool(m & 0004)
    group_read = bool(m & stat.S_IRGRP)
    return group_read


def check_files(files):
    missing_files = []
    for file in files:
        short_name = file.split("/")[-1]
        if os.path.exists(file):
            readable = group_readable(file)
            if readable:
                logger.info("%s: OK" % short_name)
            else:
                missing_files.append(file)
    if missing_files:
        logger.error("ERROR: The following files missing")
        for file in missing_files:
            logger.info("Missing %s" % file)
        sys.exit()


def delete_files_by_extension(extension):
    files = [ f for f in os.listdir(".") if f.endswith(tuple(extension))]
    for f in files:
        os.remove(f)


def delete_files(files):
    for f in files:
        os.remove(f)


def qsub_scripts(scripts):
    """ qsub scripts """
    wkdir = os.getcwd()
    for script in scripts:
        p = subprocess.Popen('ssh m0001 \"cd %s;  qsub %s\"' %
                             (wkdir, script),  shell=True,
                             stdout=subprocess.PIPE)
        output,  err = p.communicate()



def detect_cluster_job_status(completeion_file_list):
    completed = False
    for file in completeion_file_list:
        if (os.path.exists(file)):
            completed = True
        else:
            completed = False
            break
    return completed


def detect_cluster_jobs(complete_stamps):
    """ detect if job on cluster finised """
    job_status = False
    logger.info("Waiting for cluster jobs to finish!\n")
    while (not job_status):
        time.sleep(10)
        job_status = detect_cluster_job_status(complete_stamps)
    logger.info("All cluster jobs finished? %s\n" % job_status)


# def make_rscripts(patient_files, template_dir, wkdir, snv_input, cnv_input):
#     """ genearate pileup scripts, which also parses pileup results! """
#     pileup_scripts = []
#     pileup_stamps = []
#     for patient in patient_files:
#         for status in patient_files[patient]:
#             snv_positions = ".".join([patient, "vcf.snp.pos"])
#             DNA_bam = patient_files[patient][status]['DNA_bam']
#             RNA_bam = patient_files[patient][status]['RNA_bam']
#             if (DNA_bam != 'NA'):
#                 script = populate_template(patient, status, 'DNA',
#                                      snv_positions, DNA_bam, template_dir)
#                 pileup_scripts.append(script)
#                 pileup_stamp = ".".join([patient, status, 'DNA',
#                                          "pileup.complete"])
#                 pileup_stamps.append( pileup_stamp )
#             if (RNA_bam != 'NA'):
#                 script = populate_template(patient, status, 'RNA',
#                                      snv_positions, RNA_bam, template_dir)
#                 pileup_scripts.append(script)
#                 pileup_stamp = ".".join([patient, status, 'RNA',
#                                          "pileup.complete"])
#                 pileup_stamps.append( pileup_stamp )

#     return [pileup_scripts, pileup_stamps]


def populate_template(patient_status, template_dir, wkdir,
                      snv_input, cnv_input, rscript):
    rscript = ".".join([patient_status, "r"])
    jinja2_env = Environment(loader=FileSystemLoader([template_dir]),
                             trim_blocks=True)
    template = jinja2_env.get_template('expands_template.r')
    with open(rscript, 'wb') as opf:
        content = template.render(wkdir=wkdir,
                                  patient_status=patient_status,
                                  snv_input=snv_input,
                                  cnv_input=cnv_input,
                                  rscript=rscript)
        opf.write(content)
        logger.info('templated {0}'.format(rscript))
    return rscript



# def make_expands_input_file(filtered_summary, patient_files_wn):
def make_variant_dict(filtered_summary):
    """
    d = gene_variant_patients dictionary 
    tc_d = tumor content dictionary
    """
    variant_dict = dict()
    tc_dict = dict()
    with open (filtered_summary, 'r') as fh:
        records = csv.DictReader(fh,  delimiter='\t')
        for line in records:
            gene = line['gene']
            chromosome = line["chromosome"]
            pos = line["position"]
            ref = line["ref_base"]
            alt = line["alt_base"]
            DNA_n_altC = int(line["n_DNA_AltC"])
            DNA_n_af = float(line["n_DNA_AF"])
            DNA_t_refC = int(line["t_DNA_RefC"])
            DNA_t_altC = int(line["t_DNA_AltC"])
            DNA_t_af = float(line["t_DNA_AF"])
            patient = line["patient_ID"]
            # only account autosomes, and tumour af has to be > n_af
            if chromosome.isdigit() and (DNA_n_af < DNA_t_af):
                if (DNA_n_af < 0.03 or DNA_n_altC < 2):
                    alt_ploidy = 0 # somatic
                else:
                    alt_ploidy = 1 # germline
                variant = ":".join([chromosome, pos,
                                    str(DNA_t_refC), str(DNA_t_altC),
                                    str(DNA_t_af), str(alt_ploidy), gene])
                try:
                    variant_dict[patient].append(variant)
                except KeyError:
                    #print "key error!"
                    variant_dict[patient] = [variant]
    return variant_dict


def make_snv_input_file(variant_dict):
    print "Making expands tsv input files!"
    snv_header = ["chr", "startpos", "T_refC", "T_altC",
                  "AF_Tumor", "PN_B", "gene"]
    # header = ["mutation_id", "ref_counts", "var_counts","normal_cn",
    # "minor_cn", "major_cn", "variant_freq"]
    for patient in variant_dict:
        snv_file = "".join([patient, ".expands.snv"])
        with open(snv_file,  'wb') as fh:
            writer = csv.writer(fh, delimiter='\t')
            writer.writerow(snv_header)
            mutations = list(set(variant_dict[patient]))
            for mutation in mutations:
                # sl = mutation.split(":")
                writer.writerow(mutation.split(":"))


def make_expands_Rscript(patient_files_wn, template_dir):
    print "Generating expands R scripts!"
    rscripts = []
    for patient in patient_files_wn:
        for status in patient_files_wn[patient]:
            if ("normal" not in status.lower()):
                wkdir = os.getcwd() + "/"
                patient_status = "_".join([patient, status])
                directory = "".join([wkdir, patient_status])
                cnv_input = "".join([directory, ".expands.cn"])
                snv_input = "".join([directory, ".expands.snv"])
                rscript = ".".join([patient_status, "r"])
                rscripts.append(rscript)
                populate_template(patient_status, template_dir,
                                  wkdir, snv_input, cnv_input, rscript)
    return rscripts

def make_cnv_input_file(patient_files_wn):
    cn_header = ["chr", "startpos", "endpos", "CN_Estimate"]
    wkdir = os.getcwd() + "/"
    for patient in patient_files_wn:
        for status in patient_files_wn[patient]:
            if ("normal" not in status.lower()):
                blah = "_".join([patient, status])
                cn_file = "".join([wkdir, blah, ".expands.cn"])
                # copy copy number file and add header,
                # copy number needs to be absolute value
                cnv = patient_files_wn[patient][status]["cnv"]
                shutil.copyfile(cnv, cn_file)
                for line in fileinput.input(cn_file, inplace=True):
                    if fileinput.isfirstline():
                        print '\t'.join(cn_header)
                    print line,


def group_patients(patient_files):
    patient_files_wn = dict()
    patient_files_non = dict()
    for patient in patient_files:
        tmp_dict = patient_files[patient]
        if ("normal" in patient_files[patient]):
            patient_files_wn[patient] = tmp_dict
        else:
            patient_files_non[patient] = tmp_dict
    return [patient_files_wn, patient_files_non]


def run_expands(R_path, r_script):
    logger.info('running expands script: %s' % r_script)
    p = subprocess.Popen('%s/Rscript %s' % (R_path, r_script),
                         shell=True, stdout=subprocess.PIPE)
    output,  err = p.communicate()
    print output


def __main__():
    parser = argparse.ArgumentParser(description='run expands to show clonal shift in tumours!')
    parser.add_argument('-i1', '--meta_file', required=True,
                        help='specify the file, which tells vcf, bam, cnv paths')
    parser.add_argument('-i2', '--filtered_summary_file', required=True,
                        help='specify file, which lists all somatic and germline variants')
    args = vars(parser.parse_args())
    logger.debug('some debug message')
    logger.info('some info message')
    logger.warning('some warning message')
    logger.error('some error message')
    logger.critical('some critical message')
    logger.info("expands scripts starts at: %s\n" % datetime.datetime.now())
    meta_file = args['meta_file']
    filtered_summary = args['filtered_summary_file']
    template_dir = '/home/szong/projects/development/expands/'
    R_path = '/gsc/software/linux-x86_64-centos6/R-3.1.1/bin/'
    patient_files = get_files(meta_file)
    # out_dict = group_patients(patient_files)
    patient_files_wn = group_patients(patient_files)[0]
    variant_dict = make_variant_dict(filtered_summary)
    make_snv_input_file(variant_dict)
    make_cnv_input_file(patient_files_wn)
    rscripts = make_expands_Rscript(patient_files_wn, template_dir)
    for rscript in rscripts:
        print rscript
        run_expands(R_path, rscript)
    logger.info("run_expands scripts finished on: %s" %
                datetime.datetime.now())


if __name__ == '__main__':
    __main__()

