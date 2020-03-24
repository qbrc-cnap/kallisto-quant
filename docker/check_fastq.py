import subprocess as sp
import argparse
import sys
import os

R1 = 'r1'
R2 = 'r2'

def run_cmd(cmd, return_stderr=False):
    '''
    Runs a command through the shell
    '''
    p = sp.Popen(cmd, shell=True, stderr=sp.PIPE, stdout=sp.PIPE)
    stdout, stderr = p.communicate()
    if return_stderr:
        return (p.returncode, stderr.decode('utf-8'))
    return (p.returncode, stdout.decode('utf-8'))


def check_fastq_format(f):
    '''
    Runs the fastQValidator on the fastq file
    IF the file is invalid, the return code is 1 and
    the error goes to stdout.  If OK, then return code is zero.
    '''
    cmd = 'fastQValidator --file %s' % f
    rc, stdout_string = run_cmd(cmd)
    if rc == 1:
        return [stdout_string]
    return []

def check_gzip_format(f):
    '''
    gzip -t <file> has return code zero if OK
    if not, returncode is 1 and error is printed to stderr
    '''
    cmd = 'gzip -t %s' % f
    rc, stderr_string = run_cmd(cmd, return_stderr=True)
    if rc == 1:
        return [stderr_string]
    return []


def catch_very_long_reads(f, N=100, L=300):
    '''
    In case we get non-illumina reads, they will not exceed some threshold (e.g. 300bp)
    '''
    err_list = []
    zcat_cmd = 'zcat %s | head -%d' % (f, 4*N)
    rc, stdout = run_cmd(zcat_cmd)
    lines = stdout.split('\n')
        
    # iterate through the sampled sequences.  
    # We don't want to dump a ton of long sequences, so if we encounter
    # ANY in our sample, save an error message and exit the loop.
    # Thus, at most one error per fastq.
    i = 1
    while i < len(lines):
        if len(lines[i]) > L:
            return ['Fastq file (%s) had a read of length %d, '
                'which is too long for a typical Illumina read.  Failing file.' % (f, len(lines[i]))]
        i += 4
    return []


def sample_read_ids(f, N=100):

    zcat_cmd = 'zcat %s | head -%d' % (f, 4*N)
    rc, stdout = run_cmd(zcat_cmd)
    lines = stdout.split('\n')
        
    # iterate through the sampled sequences, grab the seq ID
    seq_ids = []
    for i in range(0, len(lines), 4):
        seq_ids.append(lines[i].split(' ')[0])
    return seq_ids    


def get_commandline_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r1', required=True, dest=R1)
    parser.add_argument('-r2', required=True, dest=R2)
    args = parser.parse_args()
    return vars(args)


if __name__ == '__main__':
    arg_dict = get_commandline_args()
    fastqs = [arg_dict[R1], arg_dict[R2]]

    # collect the error strings into a list, which we will eventually dump to stderr
    err_list = []

    for fastq_filepath in fastqs:
        # check that fastq in gzip:
        err_list.extend(check_gzip_format(fastq_filepath))

        # check the fastq format
        err_list.extend(check_fastq_format(fastq_filepath))

        # check that read lengths are consistent with Illumina:
        err_list.extend(catch_very_long_reads(fastq_filepath))

    # check if sorted the same way
    r1_ids = sample_read_ids(arg_dict[R1])
    r2_ids = sample_read_ids(arg_dict[R2])

    # if the format that ends with /1 or /2 to denote pairs, strip that off
    if r1_ids[0].endswith('/1'):
        r1_ids = [x[:-2] for x in r1_ids]
        r2_ids = [x[:-2] for x in r2_ids]

    # now compare the lists to see that they are equivalent
    if not r1_ids == r2_ids:
        err_list.append('The sequence identifiers suggested that the paired fastq are not properly sorted as pairs.')

    if len(err_list) > 0:
        sys.stderr.write('#####'.join(err_list)) # the 5-hash delimiter since some stderr messages can be multiline
        sys.exit(1) # need this to trigger Cromwell to fail
        
