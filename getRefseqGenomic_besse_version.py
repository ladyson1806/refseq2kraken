#!/usr/bin/env python3
"""
NAME: getKrakenFasta.py
=========

DESCRIPTION
===========

INSTALLATION
============

USAGE
=====

VERSION HISTORY
===============

0.0.1   2016/12/06    Initial version.

LICENCE
=======
2016, copyright Sebastian Schmeier (s.schmeier@gmail.com), sschmeier.com

template version: 1.6 (2016/11/09)

"""
from timeit import default_timer as timer
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np 
import sys
import os
import os.path
import argparse
import csv
import gzip
import bz2
import zipfile
import urllib
import hashlib
import multiprocessing
import subprocess
import signal
import time


__version__ = '0.0.2'
__date__ = '2022/11/16'
__email__ = 'savandara.besse@cri-paris.org'
__author__ = 'Savandara Besse'


def parse_cmdline():
    """ Parse command-line args. """
    description = 'Download fasta-genomic sequences from ncbi using rsync.'
    version = 'version %s, date %s' % (__version__, __date__)
    epilog = 'Copyright %s (%s)' % (__author__, __email__)

    parser = argparse.ArgumentParser(description=description, epilog=epilog)

    parser.add_argument('--version',
                        action='version',
                        version='%s' % (version))

    parser.add_argument('-b',
        '--branch',
        dest='str_branch',
        metavar='BRANCH',
        type=str,
        default="bacteria,viral,fungi,protozoa,archaea",
        help='Branches of organisms to download separated by comma, e.g. bacteria,fungi, etc. [default=" bacteria,viral,fungi,protozoa,archaea"]')
    
    parser.add_argument('-l',
        '--level',
        dest='str_level',
        metavar='LEVEL',
        type=str,
        default="Complete Genome",
        help='Assembly - level of genomic sequences to include, separated by comma. For example: Chromosome,Contig,Scaffold. [default="Complete Genome"]')

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
        'Number of sub-processes (concurrent downloads) to use.'+\
        ' It is only logical to not give more processes'+\
        ' than cpus/cores are available. [default: 1]')

    # if no arguments supplied print help
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args()
    return args, parser


def init_worker(tqdm_lock=None):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    if tqdm_lock is not None:
        tqdm.set_lock(tqdm_lock)
    

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


def my_func(args):
    """
    THIS IS THE ACTUAL WORKFUNCTION THAT HAS TO BE EXECUTED MULTPLE TIMES.
    This function could be distributed to the cores requested.
    # do stuff
    Here we download a file, get a status and adjust fasta-header
    return (args, res)
    """
    fname = args[0]
    dnlurl = args[1]
    dest_dir = args[2]
    rsync_cmd = f"rsync --omit-dir-times --rsync-path='/usr/bin/sudo /bin/rsync' --copy-links --partial {dnlurl} {dest_dir}"
    retcode = subprocess.call(rsync_cmd, shell=True)
    return (args, retcode)


def parse_assemblyfile(branch, genomictypes=["Complete Genome"], dest_dir='/media/UTELifeNAS/homes/Savvy/Genomes'):
    fname = f'genomes/refseq/{branch}/assembly_summary.txt'
    url = f'rsync://ftp.ncbi.nlm.nih.gov/{fname}'

    rsync_cmd = f"rsync --omit-dir-times --rsync-path='/usr/bin/sudo /bin/rsync'  --copy-links --partial {url} {dest_dir}"
    retcode = subprocess.call(rsync_cmd, shell=True)
   
    jobs = []
    # read file, extract ftp paths and download each file
    oR = csv.reader(load_file(fname), delimiter = '\t')
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
            ftp_path = a[19]
            name     = os.path.basename(ftp_path) + '_genomic.fna.gz'
            dnlurl   = os.path.join(ftp_path, name)
            #### Replace 'ftp://' by 'https://' 
            dnlurl = dnlurl.replace('https://', 'rsync://')
            jobs.append((name, dnlurl, f'/media/UTELifeNAS/homes/Savvy/Genomes/refseq/{branch}', branch))
    return jobs, retcode, d


def main():
    """ The main function. """
    args, parser = parse_cmdline()

    process_number = args.process_number
    if process_number < 1:
        parser.error('-p has to be > 0: EXIT.')

    branches = [s.strip() for s in args.str_branch.split(',')]
    types = [s.strip() for s in args.str_level.split(',')]

    job_list = []
    for branch in branches:
        job_list_branch, retcode, dStats = parse_assemblyfile(branch, types, '/media/UTELifeNAS/homes/Savvy/Genomes/')
        job_list += job_list_branch
        if args.assemblystats:
            status = dStats.keys()
            #### Change status sort function (modified by Savandara Besse)
            status = sorted(status)
            sys.stdout.write('Branch: %s\n'%branch)
            for s in status:
                sys.stdout.write('%s\t%i\n' %(s,dStats[s]))

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

    #-------------------------------------------------------------------------
    # MULTITHREADING V1 (by Sebastian Schmeier)
    #-------------------------------------------------------------------------
    # start_time = timer()  # very crude timing
    # # create pool of workers ---------------------
    # pool = Pool(processes = process_number)
    # result_list = pool.map_async(my_func, job_list, chunksize=1)  # chunksize=1 for correct progressbar
    # pool.close()  # No more work

    # jobs_total = len(job_list)
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
    #     time.sleep(2)  # wait a bit: here we test every 2 sec
    # # Finish the progress bar
    # bar_str = '=' * progress_bar_length
    # sys.stderr.write("JOBS (%s): [%s] (%i) 100%%\n" % ('0'.rjust(len(str(jobs_total))),
    #                                                    bar_str,
    #                                                    jobs_total))
    # #result_list = result_list.get()
    # end_time = timer()
    # sys.stderr.write('PROCESS-TIME: %.1f sec\nDONE.\n\n' % (end_time - start_time))
    #-------------------------------------------------------------------------

if __name__ == '__main__':
    sys.exit(main())

