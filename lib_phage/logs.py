### functions for step validations and logs writing ###

import os
import glob

def check_input_repr_prot_selection():

    # check input files integrity & if this step was already executed:
    # if it was warn about data overwrite

    """
    input.file.ok <- file.exists(input.fasta) & file.size(input.fasta)>100
      if(!input.file.ok){
        cat("Something wrong with the input file! Exiting...\n")
      }

      else{
        if(!dir.exists(output.dir)) dir.create(output.dir, recursive = T, showWarnings = F)
        mmseqs.tmp <- tempdir()
        if(dir.exists(mmseqs.tmp)) unlink(mmseqs.tmp, recursive = T); dir.create(mmseqs.tmp)
        # cat(mmseqs.tmp,"\n")
        output.filename <- sprintf("%s/%s", output.dir, output.filename)
        output.file.exists <- file.exists(output.filename) & file.size(output.filename)>100
        if(output.file.exists){
          if(overwrite){
            unlink(output.filename)
          } else{
            stop("Output file already exists! Exiting...\n")
          }
        }
    """

    return True

def validate_output_repr_prot_selection(work_dir, output_dirpath,
                                    cluster_params_min_seqid,
                                    cluster_params_sensitivity,
                                    cluster_params_coverage):

    # validate clustering and selection process - check files integrity, write to log:
    flog = open(work_dir + 'log/repr-seq.log', 'w')
    # check files present and not empty
    for f in ['/clustering.tsv', '/name-table.txt', '/repr-seqs-lengths.txt', '/repr-seqs.fa']:
        if os.path.isfile(output_dirpath + f):
            if os.path.getsize(output_dirpath + f) != 0:
                pass
            else:
                print(f + ' file empty')
                return False
        else:
            print(f + ' does not exist')
            return False

    # write validation result
    flog.write('STATUS:OK\n')

    # write params used
    flog.write('PARAMS USED:\n')
    flog.write('cluster_params_min_seqid=' + str(cluster_params_min_seqid) + '\n')
    flog.write('cluster_params_sensitivity=' + str(cluster_params_sensitivity) + '\n')
    flog.write('cluster_params_coverage=' + str(cluster_params_coverage) + '\n')

    # paths to inputs ?

    flog.close()

    print('Validation success, log file stored.')

def check_input_all_vs_all_HMM(work_dir, force=False):

    if os.path.isfile(work_dir + 'log/repr-seq.log'):
        flog = open(work_dir + 'log/repr-seq.log', 'r')
        # read status line
        status = flog.readline().split(':')[1]
        if status.strip() == 'OK':
            pass
        else:
            print('Log file reads: previous step failed. Aborting.')
    else:
        print('Failed to fetch log file, aborting.')
        return False

    # check if data already exist
    output_hhblits_dirpath = work_dir + 'intermediate/prot-families/all-by-all'

    if force == True:
        confirm = input('This will overwrite all data from this step. Proceed? [y/n]')
        if confirm == 'y':
            #!TODO clear dir
            # clear dir
            print('Clearing ' + output_hhblits_dirpath + '...')
            return True
        else:
            return False
    else:
        if len(os.listdir(output_hhblits_dirpath) ) != 0:
            print('This step was already executed. Run validation fucntion with force=True to overwrite.')
            return False

    return True

def save_params_hhblits(work_dir, n, mact, p, qid, cov):

    # validate hhblits - save params to log:
    flog = open(work_dir + 'log/hhblits.log', 'w')

    # write validation result
    flog.write('STATUS:RUNNING\n')

    # write params used
    flog.write('PARAMS USED:\n')
    flog.write('number of iterations[n]=' + str(n) + '\n')
    flog.write('posterior prob threshold for MAC realignment[mact]=' + str(mact) + '\n')
    flog.write('minimum probability in summary and alignment list[p]=' + str(p) + '\n')
    flog.write('minimum sequence identity with master sequence[qid]=' + str(qid) + '\n')
    flog.write('minimum coverage with master sequence (%)[cov]=' + str(cov) + '\n')

    # paths to inputs ?

    flog.close()

    print('Parameters saved, log file stored.')

def validate_output_hhblits(work_dir):

    # check if hhblits step complete
    flog = open(work_dir + 'log/hhblits.log', 'r')
    status = flog.readline().split(':')[-1].strip()
    params = flog.read()
    flog.close()
    if status == 'RUNNING':
        # check number of .hhr files compared to input files
        ind_seqs_dirpath       = work_dir + 'tmp/all-by-all/individual-seqs'
        output_hhblits_dirpath = work_dir + 'intermediate/prot-families/profiles'

        input_seqs_n = len(glob.glob(ind_seqs_dirpath + '/*.fa'))
        output_a3m_n = len(glob.glob(output_hhblits_dirpath + '/*.a3m'))

        if input_seqs_n == output_a3m_n:
            flog = open(work_dir + 'log/hhblits.log', 'w')
            flog.write('STATUS:OK\n')
            flog.write(params)
            flog.close()
            print('HHblits step complete. Updated status.')
            return True
        else:
            print('HHblits step not yet finished. Please wait and try again later.')
            return False

    elif status == 'OK':
        return True

def validate_create_db(work_dir):

    flog   = open(work_dir + 'log/hh-db.log', 'r')
    status = flog.readline().split(':')[-1].strip()
    flog.close()

    # write validation result
    if status == 'OK':
        return True
    else:
        return False

def validate_search_all_vs_all(work_dir):

    # check if hhblits step complete
    flog = open(work_dir + 'log/hhblits-all-vs-all.log', 'r')
    status = flog.readline().split(':')[-1].strip()
    params = flog.read()
    flog.close()
    if status == 'RUNNING':
        # check number of .hhr files compared to input files
        ind_profile_dir        = work_dir + 'intermediate/prot-families/profiles'
        output_hhblits_dirpath = work_dir + 'intermediate/prot-families/all-by-all'

        input_prof_n = len(glob.glob(ind_profile_dir + '/*.a3m'))
        output_hhr_n = len(glob.glob(output_hhblits_dirpath + '/*.hhr'))

        # first check filesize of all hhr files
        #!TODO

        if input_prof_n == output_hhr_n:
            flog = open(work_dir + 'log/hhblits-all-vs-all.log', 'w')
            flog.write('STATUS:OK\n')
            flog.write(params)
            flog.close()
            print('HHblits all-vs-all step complete. Updated status.')
            return True
        else:
            print('HHblits all-vs-all step not yet finished. Please wait and try again later.')
            return False

    elif status == 'OK':
        return True
