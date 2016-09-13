#! /usr/bin/env python

"""
This script does the following:
-> make expands input files
-> run expands
"""


import os.path, time, datetime, subprocess
import sys, glob, argparse, csv, shutil, fileinput


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


def make_variant_dict(filtered_summary):
    variant_dict = dict()
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
                    alt_ploidy = 0  # somatic
                else:
                    alt_ploidy = 1  # germline
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


def parse_args():
    parser = argparse.ArgumentParser(
        description='Run EXPANDS to detect clonal shift in tumours!')
    parser.add_argument(
        '-i1', '--meta_file', required=True,
        help='Specify the fill containing vcf, bam, cnv paths!')
    parser.add_argument(
        '-i2', '--filtered_summary_file', required=True,
        help='Specify file containing all somatic variants')
    args = parser.parse_args()
    return args


def __main__():
    logger.info("expands scripts starts at: %s\n" % datetime.datetime.now())
    args = parse_args()
    meta_file = args.meta_file
    filtered_summary = args.filtered_summary_file
    template_dir = '/home/szong/projects/development/cfanalyzer/'
    R_path = '/gsc/software/linux-x86_64-centos6/R-3.1.1/bin/'
    patient_files = get_files(meta_file)
    patient_files_wn = group_patients(patient_files)[0]
    variant_dict = make_variant_dict(filtered_summary)
    make_snv_input_file(variant_dict)
    make_cnv_input_file(patient_files_wn)
    rscripts = make_expands_Rscript(patient_files_wn, template_dir)
    # run expands jobs 
    for rscript in rscripts:
        print rscript
        run_expands(R_path, rscript)
    logger.info("run_expands scripts finished on: %s" %
                datetime.datetime.now())


if __name__ == '__main__':
    __main__()
