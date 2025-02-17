#!/usr/bin/env python3
"""
NAME: getKrakenFna.py
=========

DESCRIPTION
===========

INSTALLATION
============

USAGE
=====

VERSION HISTORY
===============

0.0.2   2017/01/20    Allowed for gziped output files
0.0.1   2016/12/06    Initial version.

LICENCE
=======
2016, copyright Sebastian Schmeier (s.schmeier@gmail.com), sschmeier.com

template version: 1.6 (2016/11/09)
"""
from timeit import default_timer as timer
from multiprocessing import Pool
from tqdm import tqdm 
from Bio import SeqIO # non-standard lib BioPython
import numpy as np 
import sys
import os
import os.path
import argparse
import csv
import gzip
import bz2
import multiprocessing
import signal
import subprocess
import zipfile
import urllib
import hashlib
import time


__version__ = '0.0.3'
__date__ = '2022/11/16'
__email__ = 'savandara.besse@cri-paris.org'
__author__ = 'Savandara Besse'

def init_worker(tqdm_lock=None):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    if tqdm_lock is not None:
        tqdm.set_lock(tqdm_lock)
    

def parse_cmdline():
    """ Parse command-line args. """
    ## parse cmd-line -----------------------------------------------------------
    description = "Process fasta-genomic sequences from NCBI-refseq for inclusion in a KrakenDB. A branch's assembly_summary.txt as well as all of its *_genomic.fna.gz files are assumed to be in a sub-directory called like the branch, e.g. ./genomes/refseq/bacteria"
    version = 'version %s, date %s' % (__version__, __date__)
    epilog = 'Copyright %s (%s)' % (__author__, __email__)

    parser = argparse.ArgumentParser(description=description, epilog=epilog)

    parser.add_argument('--version',
                        action='version',
                        version='%s' % (version))

    parser.add_argument(
        dest='str_kraken',
        metavar='KrakenDB-DIR',
        type=str,
        help='Directory name for the processed fasta-files for kraken. Will be created if not found.')

    parser.add_argument('-b',
        '--branch',
        dest='str_branch',
        metavar='BRANCH',
        type=str,
        default="bacteria,viral,fungi,protozoa,archaea",
        help='Branches of organisms to download separated by comma, e.g. bacteria, fungi, etc. [default=" bacteria,viral,fungi,protozoa,archaea"]')
    
    parser.add_argument('-l',
        '--level',
        dest='str_level',
        metavar='LEVEL',
        type=str,
        default="Complete Genome",
        help='Assembly - level of genomic sequences to include, separated by comma. For example: Chromosome,Contig,Scaffold. [default="Complete Genome"]')
    parser.add_argument('-d',
        '--dir',
        dest='str_dir',
        metavar='DIRECTORY',
        type=str,
        default="/media/UTELifeNAS/homes/Savvy/Genomes/refseq/",
        help='Base directory for refseq fasta-files. Here, we assume sub-directories for branches, e.g. bacteria etc. [default="/media/UTELifeNAS/homes/Savvy/Genomes/refseq/"]')

    parser.add_argument('-z',
        '--gzip',
        dest='gzip',
        default=False,
        action='store_true',
        help='Create gzip-ed output files. (process takes a lot longer)')

    parser.add_argument('-a',
        '--assembly',
        dest='assemblystats',
        default=False,
        action='store_true',
        help='Print assembly stats for branches and exits.')

    group1 = parser.add_argument_group('Threading',
                                       'Multithreading arguments:')

    group1.add_argument(
        '-p',
        '--processes',
        metavar='INT',
        type=int,
        dest='process_number',
        default=1,
        help=
        'Number of concurrent sub-processes to use.'+\
        ' It is only logical to not give more processes'+\
        ' than cpus/cores are available. [default: 1]')

    # if no arguments supplied print help
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args()
    return args, parser

def load_file(filename):
    """ LOADING FILES """
    if filename in ['-', 'stdin']:
        filehandle = sys.stdin
    elif filename.split('.')[-1] == 'gz':
        filehandle = gzip.open(filename)
    elif filename.split('.')[-1] == 'bz2':
        filehandle = bz2.BZFile(filename)
    elif filename.split('.')[-1] == 'zip':
        filehandle = zipfile.Zipfile(filename)
    else:
        filehandle = open(filename)
    return filehandle

def unzip_file(filename):
    if filename.split('.')[-1] == 'gz':
        gunzip_cmd = f'gunzip {filename}' 
        subprocess.run(gunzip_cmd, shell=True)
        filehandle = filename.replace('.gz', '')
    else:
        filehandle = filename
    return filehandle


def my_func(args):
    """
    THIS IS THE ACCTUAL WORK FUNCTION THAT HAS TO BE EXECUTED MULTPLE TIMES.
    This function could be distributed to the cores requested.
    # do stuff
    args = (taxid, infile-path, outfile-path)
    
    return (args, res)
    """
    taxid = args[0]
    infilepath = args[1]
    outfilename = args[2]
    if not os.path.isfile(infilepath):
        sys.stderr.write('%s not found. SKIP\n'%(infilepath))
        return (args, 0)
    else:
        unzip_file(infilepath)
    
    return (args, 1)


def parse_assemblyfile(branch, genomictypes=["Complete Genome"], dirpath='/media/UTELifeNAS/homes/Savvy/Genomes', krakendir='/media/UTELifeNAS/homes/Savvy/BBDuk', gzip=False):
    basedir = os.path.join(dirpath, branch)
    fname = 'assembly_summary.txt'
    krakendir = os.path.join(krakendir, branch)

    if not os.path.isfile(os.path.join(basedir,fname)):
        sys.stderr.write("ERROR: '%s' not for branch '%s' at %s\nUse 'python getRefseqGenomic.py -b %s -p 4' to download fasta-sequences and assembly-summary.txt.\nEXIT\n\n"
                             % (fname, branch, os.path.join(basedir,fname), branch))
        sys.exit(1)
    else:
        jobs = []
        # read file, extract ftp paths and download each file
        oR = csv.reader(load_file(os.path.join(basedir,fname)), delimiter = '\t')
        d = {}
        for a in oR:
            try:
                assembly_level =  a[11]
            except:
                continue

            version_status =  a[10]
            if version_status != 'latest':
                continue
        
            if assembly_level == 'assembly_level':
                continue
            
            d[assembly_level] = d.get(assembly_level, 0) + 1

            if assembly_level in genomictypes:
                name     = os.path.basename(a[19]) + '_genomic.fna.gz'
                filepath = os.path.join(basedir, name)
                taxid    = a[5]

                if gzip:
                    fnameTax = name.replace('.fna.gz', '.tax.fna.gz')  # store gziped files
                else:
                    fnameTax = name.replace('.fna.gz', '.tax.fna')

                jobs.append((taxid, filepath, os.path.join(krakendir, fnameTax)))
        
    return jobs, d


def main():
    """ The main function. """
    args, parser = parse_cmdline()

    branches = [s.strip() for s in args.str_branch.split(',')]
    types = [s.strip() for s in args.str_level.split(',')]
    dirpath = args.str_dir
    process_number = args.process_number
    if process_number < 1:
        parser.error('-p has to be > 0: EXIT.')

    job_list = []
    for branch in branches:
        job_list_br, dStats = parse_assemblyfile(branch,
                                                 types,
                                                 dirpath,
                                                 args.str_kraken,
                                                 args.gzip)
        job_list += job_list_br
        if args.assemblystats:
            status = dStats.keys()
            status = sorted(status)
            sys.stdout.write('Branch: %s\n'%branch)
            for s in status:
                sys.stdout.write('%s\t%i\n' %(s,dStats[s]))
        else:
            # prepare directory structure for results
            if not os.path.exists(os.path.join(args.str_kraken, branch)):
                sys.stdout.write(f'Build database for {branch}')

    # exit if only stats should be displayed
    if args.assemblystats:
        return

    #-------------------------------------------------------------------------
    # MULTITHREADING V2 (by Savandara Besse)
    #-------------------------------------------------------------------------
    p = multiprocessing.Pool(initializer=init_worker, initargs=(tqdm.get_lock(),), processes=process_number)
    try:
        pbar = tqdm(job_list, maxinterval=1.0, miniters=1, desc="Completed jobs: ", bar_format="{desc}:{percentage:3.0f}%|{bar}|")
        for _, result in enumerate(p.imap_unordered(my_func, job_list, chunksize=1)):
            pbar.update(1)  # Everytime the iteration finishes, update the global progress bar

        pbar.close()
        p.close()
        p.join()
    except KeyboardInterrupt:
        print("KeyboardInterrupt, terminating workers.")
        pbar.close()
        p.terminate()
        p.join()
        time.sleep(np.random.uniform(low=0.0, high=2.0, size=None))
        exit(1)

    return

    # #-------------------------------------------------------------------------
    # # MULTITHREADING V1 (by Sebastian Schmeier)
    # #-------------------------------------------------------------------------
    # # For timing
    # start_time = timer()  # very crude timing
    # # create pool of workers ---------------------
    # pool = Pool(processes=process_number)

    # # "chunksize"" usually only makes a noticable performance
    # # difference for very large iterables
    # # Here I set it to 1 to get the progress bar working nicly
    # # Otherwise it will not give the correct number of processes left
    # # to process but rather the chunksize number.
    # chunksize = 1
    # result_list = pool.map_async(my_func, job_list, chunksize=chunksize)
    # pool.close()  # No more work

    # jobs_total = len(job_list)
    # # Progress bar
    # #==============================
    # # This can be changed to make progressbar bigger or smaller
    # progress_bar_length = 50
    # #==============================
    # while not result_list.ready():
    #     num_not_done = result_list._number_left
    #     num_done = jobs_total - num_not_done
    #     num_bar_done = num_done * progress_bar_length / jobs_total
    #     bar_str = ('=' * num_bar_done).ljust(progress_bar_length)
    #     percent = int(num_done * 100 / jobs_total)
    #     sys.stderr.write("JOBS (%s): [%s] (%s) %s%%\r" % (str(num_not_done).rjust(len(str(jobs_total))),
    #                                                       bar_str,
    #                                                       str(num_done).rjust(len(str(jobs_total))),
    #                                                       str(percent).rjust(3)))

    #     sys.stderr.flush()
    #     time.sleep(1)  # wait a bit: here we test every sec
    # # Finish the progress bar
    # bar_str = '=' * progress_bar_length
    # sys.stderr.write("JOBS (%s): [%s] (%i) 100%%\n" % ('0'.rjust(len(str(jobs_total))),
    #                                                    bar_str,
    #                                                    jobs_total))
    # #result_list = result_list.get()
    # # --------------------------------------------

    # end_time = timer()
    # sys.stderr.write('PROCESS-TIME: %.1f sec\nDONE.\n\n' % (end_time - start_time))

def merge_fasta(dirpath, bbdukdir):
    args, parser = parse_cmdline()
    branch = args.str_branch
    base_folder = os.path.join(dirpath, branch)

    fasta_list = [ os.path.join(base_folder, fasta_file) for fasta_file in os.listdir(base_folder) if '.fna' in fasta_file ] 
    for i in range(len(fasta_list)):
        fasta_file = fasta_list[i]
        if i == 0 : 
            merge_cmd = f'cat {fasta_file} > {bbdukdir}/{branch}_database.fna'
        else : 
            merge_cmd = f'cat {fasta_file} >> {bbdukdir}/{branch}_database.fna'
        subprocess.run(merge_cmd, shell=True)

        gzip_cmd = f'gzip {fasta_file} &'
        subprocess.run(gzip_cmd, shell=True)
    

    gzip_cmd = f'gzip *.fna'


if __name__ == '__main__':
    main()
    merge_fasta('/media/UTELifeNAS/homes/Savvy/Genomes/refseq/', '/media/UTELifeNAS/homes/Savvy/BBDuk/')
    sys.exit()

