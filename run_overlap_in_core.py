#!/usr/bin/env python
import argparse
import os
import re
import subprocess
import sys

def get_parser():

    parser = argparse.ArgumentParser()
    parser.add_argument('--ovl_partitions', required=True)
    parser.add_argument('--hashbits', required=True)
    parser.add_argument('--ca_bin_dir', required=True)
    parser.add_argument('--gkp_store', required=True)
    parser.add_argument('--nmers_fasta', required=True)
    parser.add_argument('--output_dir', required=True)
    parser.add_argument('--do_partial', action='store_true')
    
    return parser

def setup_env():
    """This is very sad."""

    os.environ['AS_OVL_ERROR_RATE'] = '0.06'
    os.environ['AS_CNS_ERROR_RATE'] = '0.25'
    os.environ['AS_CGW_ERROR_RATE'] = '0.25'
    os.environ['AS_OVERLAP_MIN_LEN'] = '40'
    os.environ['AS_READ_MIN_LEN'] = '500'

def run_ovl(opt_vals, ca_bin_dir, hashbits, nmers_fasta, gkp_store, output_dir, do_partial):

    bounds = re.findall("\d+", opt_vals)
    out_file = os.path.join(output_dir,
                            "ca_{j}.ovb.gz".format(j='_'.join([str(k) for k in bounds])))
    
    if do_partial:
        cmd_template = ("{ca_bin}/overlapInCore -G --hashbits {hashbits} --hashload 0.75 -t 1 {opt} "
                        "-k 14 -k {nmers_fasta} -o {out_file} {gkp}")
    else:
        cmd_template = ("{ca_bin}/overlapInCore --hashbits {hashbits} --hashload 0.75 -t 1 {opt} "
                        "-k 14 -k {nmers_fasta} -o {out_file} {gkp}")

    cmd = cmd_template.format(ca_bin=ca_bin_dir,hashbits=hashbits, opt=opt_vals,
                              nmers_fasta=nmers_fasta, out_file=out_file,
                              gkp=gkp_store)
    print "Running", cmd
    proc = subprocess.Popen(cmd, shell=True)
    stdout, stderr = proc.communicate()
    
    return proc.returncode

def main():

    parser = get_parser()
    args = parser.parse_args()
    
    partitions = [k.strip() for k in file(args.ovl_partitions).readlines()]

    setup_env()

    for partition in partitions:
        ret = run_ovl(partition, args.ca_bin_dir, args.hashbits, args.nmers_fasta,
                      args.gkp_store, args.output_dir, args.do_partial)
        if ret != 0:
            print "Overlap job failed:", ovl_input
            sys.exit(1)


if __name__ == '__main__':
    main()
