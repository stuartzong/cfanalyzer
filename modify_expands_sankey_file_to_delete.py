#! /usr/bin/env python

import os, stat, os.path, time, datetime, subprocess
import re, sys, glob, argparse, csv, shutil, fileinput
from collections import defaultdict
from pprint import pprint
from itertools import islice
import ConfigParser

def __main__():
    #print "Quality and somatic filtering script starts at: %s\n" % datetime.datetime.now()     
    parser = argparse.ArgumentParser(description='Filter variants based on qulaity and somatic filters')
    parser.add_argument('-i','--input_file', help='specify sankey input file, which needs to be modified', required=True)
    args = vars(parser.parse_args())

    # variant_input_file = "bam_vcf_cnv_path.txt" 
    old_sankey = args['input_file']
    new_sankey = ".".join([old_sankey, "modified" ])

    modify_sankey(old_sankey, new_sankey)
    #clonal_change(old_sankey)

def clonal_change(old_sankey):
    d = dict
    with open (old_sankey, 'r') as handle:
            records = csv.DictReader(handle,  delimiter='\t')
            headers = records.fieldnames
            primary_clones = []
            refractory_clones = []
            for line in records:
                primary_sp = line['primary_sp']
                mutation = line["mutation"]
                gene = mutation.split('_')[0]
                refractory_sp = line["refractory_sp"]
                primary_clones.append(primary_sp)
                refractory_clones.append(refractory_sp)
            num_priclones = len(set(primary_clones)) -1
            num_refclones = len(set(refractory_clones)) -1
            #print "number of primary and refractory clones are:" 
            #print "%s\t%s"  % (old_sankey, '->'.join([str(num_priclones), str(num_refclones)]))
            print "%s\t%s\t%s"  % (old_sankey, num_priclones, num_refclones)
            if ('MODIFIER' not in mutation):
                content = [line[i] for i in headers]
                content = [primary_sp, gene, refractory_sp, value]
                print "This is a coding mutation!"
                writer.writerow(content)
            else: #make a modifier dictionary
                key = ':'.join([primary_sp, refractory_sp])
                try:
                    d[key].append(mutation)
                except KeyError:
                    d[key] = [mutation]


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
if __name__ == '__main__':
    __main__()

