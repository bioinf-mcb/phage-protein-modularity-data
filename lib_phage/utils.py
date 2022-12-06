### various pipeline functions ###

import os
import glob
import pandas as pd
import numpy as np
from subprocess import call
from Bio import SeqIO
from csb.bio.io import HHOutputParser
from shutil import copyfile
# from csb.bio.io.hhpred import HHOutputParser

def setup_dir_tree(work_dir):

	"""Setup all directories used in the pipeline."""

	### function should be commented with dir tree description ###

	### tmp dirs
	for d in ['tmp', 'tmp/repr-proteins', 'tmp/mmseqs', 'tmp/all-by-all',
			  'tmp/all-by-all/individual-seqs', 'tmp/parse', 'tmp/prot-families',
			  'tmp/prot-families/pair_table_chunks']:
		try:
			os.mkdir(work_dir + d)
		except FileExistsError:
			print('/' + d + '/ directory already set up')

	### inputs
	for d in ['input', 'input/phanotate', 'input/coding-seqs']:
		try:
			os.mkdir(work_dir + d)
		except FileExistsError:
			print('/' + d + '/ directory already set up')

	### outputs
	for d in ['output', 'output/prot-families', 'output/prot-families/representative',
			  'output/prot-families/all-by-all', '/output/prot-families/functional',
			  '/output/prot-families/families',
			  'output/prot-families/all-by-all/hhblits', 'output/prot-families/all-by-all/mmseqs']:
		try:
			os.mkdir(work_dir + d)
		except FileExistsError:
			print('/' + d + '/ directory already set up')

	### intermediates
	for d in ['intermediate', 'intermediate/prot-families', 'intermediate/prot-families/profiles',
			  'intermediate/prot-families/profiles/hhblits', 'intermediate/prot-families/functional',
			  'intermediate/prot-families/functional/hhrs',
			  'intermediate/prot-families/profiles/mmseqs', 'intermediate/prot-families/all-by-all',
			  'intermediate/prot-families/all-by-all/hhblits',
			  'intermediate/prot-families/all-by-all/mmseqs', 'intermediate/prot-families/db',
			  'intermediate/prot-families/db/hhblits',
			  'intermediate/prot-families/db/mmseqs',
			  'intermediate/prot-families/families']:
		try:
			os.mkdir(work_dir + d)
		except FileExistsError:
			print('/' + d + '/ directory already set up')

	### logs
	try:
		os.mkdir(work_dir + 'log')
	except FileExistsError:
		print('/log/ directory already set up')

def fetch_and_rename_protein_ids(work_dir, clustering_filepath, cds_all_filepath):

	"""Rename fasta files from NCBI ids to project internal ids."""

	clustering_results = pd.read_csv(clustering_filepath, sep='\t', header=None)

	# get unique protein ids
	repr_prot_names = clustering_results[0].unique()

	# !!! minimize number of proteins for quick code check
	# repr_prot_names = repr_prot_names[:16]
	# print('DEVEL: Restricting input to 16 proteins for fast calculation')

	repr_prot_names_filepath = '{}/repr-names'.format(work_dir + 'tmp/repr-proteins')

	f = open(repr_prot_names_filepath, 'w')
	for prot in repr_prot_names:
		f.write(prot + '\n')
	f.close()

	no_repr_prot = len(repr_prot_names)


	repr_prot_seqs_filepath = '{}/repr-seqs.fa'.format(work_dir + 'tmp/repr-proteins')
	extract_repr_cmd = 'seqtk subseq {} {} > {}'.format(cds_all_filepath,
														repr_prot_names_filepath,
														repr_prot_seqs_filepath)
	call(extract_repr_cmd, shell=True)

	# create table mapping reprseq-like id to ncbi id
	name_table_filepath = '{}/name-table.txt'.format(work_dir + 'output/prot-families/representative')

	f = open(name_table_filepath, 'w')
	f.write('repr.name,prot.name\n')
	protid2reprid = {}
	for repr, prot in enumerate(repr_prot_names):
		reprid = 'reprseq{:0>{id_width}}'.format(repr+1, id_width=len(str(no_repr_prot)))
		f.write(reprid + ',' + prot + '\n')
		protid2reprid[prot] = reprid
	f.close()

	# rename proteins in fasta file and store file with seqs lengths
	fasta                         = SeqIO.parse(repr_prot_seqs_filepath, "fasta")
	repr_fasta                    = []
	repr_prot_seqs_filepath_final = '{}/repr-seqs.fa'.format(work_dir + 'output/prot-families/representative')
	repr_prot_lens_filepath_final = '{}/repr-seqs-lengths.txt'.format(work_dir + 'output/prot-families/representative')
	repr_prot_lens = pd.DataFrame({'name':[], 'length':[]})

	for record in fasta:
		record.id          = protid2reprid[record.id]
		record.name        = ''
		record.description = ''
		repr_prot_lens     = repr_prot_lens.append({'name':record.id, 'length':len(record)}, ignore_index=True)
		repr_fasta.append(record)
	SeqIO.write(repr_fasta, repr_prot_seqs_filepath_final, "fasta")

	repr_prot_lens['length'] = repr_prot_lens['length'].astype('int')
	repr_prot_lens           = repr_prot_lens.sort_values('name')
	repr_prot_lens.to_csv(repr_prot_lens_filepath_final, index=False)

	return no_repr_prot, name_table_filepath

def build_hhr_table(work_dir, run_mode):

	"""Build a table of results from hhr files."""

	output_hhblits_dirpath = work_dir + 'intermediate/prot-families/all-by-all/' + run_mode

	name_table_filepath = '{}/name-table.txt'.format(work_dir + 'output/prot-families/representative')
	name_table          = pd.read_csv(name_table_filepath)
	prot_n              = len(name_table)

	hhr_table_filpath   =  '{}/table-hhr.txt'.format(work_dir + 'output/prot-families/all-by-all/' + run_mode)
	ftable              = open(hhr_table_filpath, 'w')
	ftable.write('qname,qstart,qend,qlength,sname,sstart,send,slength,pident,bitscore,eval,prob,pval\n') # write header

	for fhhr in sorted(glob.glob(output_hhblits_dirpath + '/*.hhr')):
		qname    = fhhr.split('/')[-1].split('.')[0]
		parser   = HHOutputParser()
		try:
			hit_list = parser.parse_file(fhhr)
			for hit in hit_list:
				record = ','.join([ str(i) for i in [qname, hit.qstart, hit.qend,
								   hit.qlength, hit.id, hit.start, hit.end, hit.slength,
								   int(hit.identity), hit.score, hit.evalue, (hit.probability * 100),
								   hit.pvalue]])
				ftable.write(record + '\n')
		except:
			print('Alignment error in', fhhr)
	ftable.close()

def build_hhr_table_dbs(work_dir, run_mode, db_name):

	"""Build a table of results from hhr files obtained from profiles vs external db"""

	output_hhblits_dirpath = work_dir + '/intermediate/prot-families/functional/hhrs/{}'.format(db_name)

	hhr_table_filpath   = '{}/hhblits-{}.txt'.format(work_dir + 'intermediate/prot-families/functional', db_name)
	ftable              = open(hhr_table_filpath, 'w')
	ftable.write('qname,qstart,qend,qlength,sname,sstart,send,slength,pident,bitscore,eval,prob,pval,annot\n') # write header

	for fhhr in sorted(glob.glob(output_hhblits_dirpath + '/*.hhr')):
		qname    = fhhr.split('/')[-1].split('.')[0]
		parser   = HHOutputParser()
		try:
			hit_list = parser.parse_file(fhhr)
			for hit in hit_list:
				record = ','.join([ str(i) for i in [qname, hit.qstart, hit.qend,
							   hit.qlength, hit.id, hit.start, hit.end, hit.slength,
							   int(hit.identity), hit.score, hit.evalue, (hit.probability * 100),
							   hit.pvalue, hit.name.replace(',', ' ')]])
				ftable.write(record + '\n')
		except:
			print('Alignment error in', fhhr)
	ftable.close()

def parse_hhr_single_file(protein_hhr_path, ecf=False):
	"""Parse results from single hhblits file and store it in dataframe.

	Parameters:
	-----------
	protein_hhr_path : string
		Path to the hhblits results file (.hhr).
	ecf : bool
		Switch for parsing ECFs and other files. If on, it will add 'ecf_'
		identifier to the hit id.

	Returns:
	--------
	pd.DataFrame
		Dataframe with parsed results from one file.
	"""
	parser    = HHOutputParser()
	qname     = protein_hhr_path.split('/')[-1].split('.')[0]
	hits_data = {}

	try:
		hits = parser.parse_file(protein_hhr_path)
		for n, hit in enumerate(hits):
			if ecf:
				hit_id = 'ecf_' + hit.id
			else:
				hit_id = hit.id
			record = [ str(i) for i in [qname, hit.qstart, hit.qend,
					   hit.qlength, hit_id, hit.start, hit.end, hit.slength,
					   int(hit.identity), hit.score, hit.evalue,
					   (hit.probability * 100),
					   hit.pvalue, hit.name.replace(',', ' ')] ]
			hits_data[n] = record
	except:
		print('Alignment error in', protein_hhr_path)

	protein_hhr_parsed = pd.DataFrame.from_dict(hits_data, orient='index',
									  columns=['qname', 'qstart', 'qend',
									  'qlength', 'sname', 'sstart', 'send',
									  'slength', 'pident', 'bitscore', 'eval',
									  'prob', 'pval', 'annot'])

	return protein_hhr_parsed

def create_bash_script_to_parse_hhr_results(work_dir, pipeline_env_path,
											lib_pp_path, ecf=False):
	"""Create bash script to run hhr results parser in a parallel manner.

	Parameters:
	-----------
	work_dir : string
		Path to the working directory.
	pipeline_env_path : string
		Path to virtualenv where phage pipeline code is stored. It needs to be
		communicated to the script to activate virtualenv in bash script.
	lib_pp_path : string
		Path to pipeline code. It needs to be communicated to the Python
		script so the parsing functions will be loaded.
	ecf : bool
		Switch for parsing ECFs and other files. If on, it will add 'ecf_'
		identifier to the hit id.

	Returns:
	--------
	string
		Absolute path to the created bash script.
	"""
	script_filepath =  work_dir + 'tmp/parse/helper-parse-db.sh'

	cmd = '#!/bin/bash\n\n'
	cmd += 'FILE=$(basename "${1}")\nFILE=${FILE%.*}\n'
	cmd += 'source ' + pipeline_env_path + '/bin/activate\n'
	cmd += 'python3 ' + pipeline_env_path \
	+ '/phage-pipeline/lib_phage/run_parsing_db_single_protein.py ${1} ' \
	+ str(ecf) + ' ' + lib_pp_path + '\n'

	with open(script_filepath, 'w') as file_sh:
		file_sh.write(cmd)

	return script_filepath

def run_parsing_with_bash(work_dir, n_cores, db_name, run_mode=None,
						  wildcard='reprseq*.hhr'):
	"""Run previously prepared bash script that will run parsing on all
	protein db scan results in dataset in a parallel manner. Before run clean
	directory from results of previous parsing attempts if any.

	Parameters:
	-----------
	work_dir : string
		Path to the working directory.
	n_cores : int
		Maximum number of cores to be used in parallel run.
	db_name : string
		Name of the database results to parse. If None then assume it is
		all-vs-all run and parse files from all-vs-all results dir.
	run_mode : string
		Mode of run [hhblits/mmseqs]. Needed to set correct input path. Only
		used when parsing results of all-vs-all.
	wildcard : string
		Wildcard that will be used by script to find all hhr files to parse.

	Returns:
	--------
	"""
	script_filepath = work_dir + 'tmp/parse/helper-parse-db.sh'
	if db_name:
		results_dirpath = work_dir \
		+ '/intermediate/prot-families/functional/hhrs/{}'.format(db_name)
	else:
		results_dirpath = work_dir \
		+ '/intermediate/prot-families/all-by-all/{}'.format(run_mode)

	# clean previous parsing partial results
	if len(glob.glob(results_dirpath + '*.txt')) > 0:
		print('Clearing previous parsing partial files...')
		call('rm  ' + results_dirpath + '*.txt', shell=True)
		print('Cleared.')

	# set execute mode
	call('chmod a+x ' + script_filepath, shell=True)

	# run script
	run_cmd = 'nohup find {} -name "{}" | xargs -P {} -n 1 {} \
	&'.format(results_dirpath, wildcard, n_cores, script_filepath)
	call(run_cmd, shell=True)

def concatenate_parsing_results(work_dir, db_name):
	"""Check if all results were processed with parser and if yes then
	colect all the resulting dataframes into summary files. If not all expected
	output files are present function quits.

	Parameters:
	-----------
	work_dir : string
		Path to the working directory.
	db_name : string
		Name of the database results to parse.

	Returns:
	--------
	bool
		Status of the execution. Returns False if still waiting for processing
		to complete. Returns True if process complete and files were successfuly
		concatenated.
	"""
	# check if all expected output files are present
	results_dirpath = work_dir \
	+ '/intermediate/prot-families/functional/hhrs/{}'.format(db_name)
	nb_expected_files  = len(glob.glob(results_dirpath + '/*.hhr'))
	nb_complete_files  = len(glob.glob(results_dirpath + '/*.txt'))
	hhr_table_filpath  = '{}/hhblits-{}.txt'.format(work_dir
	+ 'intermediate/prot-families/functional', db_name)

	if nb_expected_files == nb_complete_files:
		print('All files processed. Concatenating...')
		results_all_part = []
		for prot_result in glob.glob(results_dirpath + '/*.txt'):
			results_all_part.append(pd.read_csv(prot_result))
		results_all = pd.concat(results_all_part)
		results_all.to_csv(hhr_table_filpath, index=False)
		print('Done.')
		return True
	else:
		print('Not all files processed. Check again later.')
		return False

def clean_clustering_partial_data(work_dir, db_name):
	"""Clean directory with single protein dataframes. Saves space.

	Parameters:
	-----------
	work_dir : string
		Path to the working directory.

	Returns:
	--------
	"""
	results_dirpath = work_dir \
	+ '/intermediate/prot-families/functional/hhrs/{}'.format(db_name)

	if len(glob.glob(results_dirpath + '/*.txt')) > 0:
		print('Clearing parsing partial files...')
		call('rm  ' + results_dirpath + '/*.txt', shell=True)
		print('Cleared.')

def setup_paths():
	'''load existing configuration from file or create new if no file'''

	data_dir  = ''
	bin_paths = {}
	temp_dir  = ''

	# check if file config.py exists
	# if True: load paths from there to variables
	# if False: create file and ask user to state paths
	# return loaded/created data

	return data_dir, bin_paths, temp_dir

def get_ecf_hits(csv):
	df = pd.read_csv(csv)
	return df[['sname', 'qstart', 'qend']].apply(lambda row:"_".join([str(i) for i in row]), axis=1).tolist()

def clean_msa(msa):
	"""Remove columns in MSA where only query sequence has amino acid."""

	msa_clean = []

	msa_cols = [ [ aa for aa in seq ] for seq in msa ]
	msa_arr  = np.array(msa_cols)

	# take each column (alignment position) one by one and iterate over
	# ommit first row (query row)
	for pos in msa_arr.T:
		for aa in pos[1:]:
			if aa != '-':
				msa_clean.append(pos)
				break

	# transpose beack to get original array orientation
	msa_clean = np.array(msa_clean)
	msa_clean = msa_clean.T

	return msa_clean

def generate_msa(hhr_file, queryseq, hitslist, maxevalue=1e-3, ident_cut=0.5, qcov_cut=0.5, eval_cut=1e-3):

	assert len(hitslist) > 0, 'provide at least one hit id in `hitlist`'

	fasta     = [queryseq]
	hit_order = []

	for hit in HHOutputParser(alignments=True).parse_file(hhr_file):

		hit_id = f'{hit.id}_{hit.qstart}_{hit.qend}'
		if not hit_id in hitslist: continue

		query_cov = 1.*len(hit.alignment.subject.replace('-', ''))
		if hit.identity < ident_cut: continue
		if query_cov / len(queryseq) < qcov_cut: continue
		if hit.evalue > eval_cut: continue

		hit_order.append(hit.id)
		temp  = ''
		mpos  = 0
		sbjct = "-"*(hit.qstart-1) + hit.alignment.subject

		# for each aa in sbjct
		for i in range(len(sbjct)):

			# no insertion at this position
			if i - hit.qstart + 1 < 0 or hit.alignment.query[i - hit.qstart+ 1 ] != '-':
				while fasta[0][mpos] == '-':
					mpos = mpos + 1
					temp = temp + '-'
				temp = temp + str(sbjct[i])
				mpos = mpos + 1

			# insertion present
			else:
				if fasta[0][mpos] != '-':
					for f in range(len(fasta)): # we need to add a gap
						fasta[f] = fasta[f][:mpos] + "-" + fasta[f][mpos:]
				temp = temp + str(sbjct[i])
				mpos = mpos + 1

		fasta.append(temp)

	# fill gaps at the N terminus
	for f in range(len(fasta)):
		if len(fasta[f])<len(fasta[0]):
			fasta[f]=fasta[f]+ "-" * (len(fasta[0])-len(fasta[f]))

	return fasta, hit_order

def process_phanotate_output(phanotate_filepath, work_dir):
	"""Process raw output from Phanotate into compressed AA fasta file that is used further in pipeline."""

	# copy file into input dir
	phanotate_filepath_work_dir = work_dir + 'input/phanotate/cds-nt.fa'
	cds_aa_filepath             = work_dir + 'input/coding-seqs/cds-aa.fa'
	copyfile(phanotate_filepath, phanotate_filepath_work_dir)

	# translate
	fasta_nt = SeqIO.parse(phanotate_filepath_work_dir, "fasta")
	fasta_aa = []

	# iterate sequences and translate each
	for record in fasta_nt:
		record.seq         = record.seq.translate()[:-1] # remove stop codon asterix
		record.name        = '' # clean header from Phanotate data
		record.description = '' # clean header from Phanotate data
		fasta_aa.append(record)
	SeqIO.write(fasta_aa, cds_aa_filepath, "fasta")

	# compress
	call('bgzip {}'.format(cds_aa_filepath), shell=True)

	print('All fasta translated. File compressed.')

def create_reprseq_profile_from_clustering(clustering_filepath, clustering_msa_filepath,
										   cds_all_filepath, profile_outdir):

	# load clusters and check one by one if bigger than 1
	# if yes: get msa for cluster and save it as a3m
	# if no: get sequence for singleton cluster and save as a3m
	# needs mapping of ncbi to reprseqs (how to handle clusteres?)

	# open clustering file and read clusters
	fclust = open(clustering_filepath, 'r')

	clusters = {}
	for record in fclust:
		cluster = record.split('\t')[0]
		node    = record.split('\t')[1].strip()
		if cluster in clusters.keys():
			clusters[cluster].append(node)
		else:
			clusters[cluster] = [node]

	fclust.close()

	# print(clustering_msa_filepath)

	# clean MSA file (raw file has some extra white chars that interfere with Biopython)
	clust_msa_clean = open(clustering_msa_filepath + '_clean', 'w')
	clust_msa_raw = open(clustering_msa_filepath, 'r')
	for line in clust_msa_raw:
		if '>' in line:
			clust_msa_clean.write('>' + line.split('>')[1])
		else:
			clust_msa_clean.write(line)
	clust_msa_raw.close()
	clust_msa_clean.close()

	# parse MSA file
	clust_msa_seqs = {}
	clust_msa      = SeqIO.parse(clustering_msa_filepath + '_clean', "fasta")
	for record in clust_msa:
		clust_msa_seqs[record.id] = record.seq

	# parse all seqs file
	cds_all = {}
	cds     = SeqIO.parse(cds_all_filepath, "fasta")
	for record in cds:
		cds_all[record.id] = record.seq

	# save each MSA as an a3m profile
	nb_clusters = len(clusters.items())
	for repr, (cluster, nodes) in enumerate(clusters.items()):
		reprid = 'reprseq{:0>{id_width}}'.format(repr+1, id_width=len(str(nb_clusters)))
		with open(profile_outdir + '/' + reprid + '.a3m', 'w') as fprofile:
			if len(nodes) > 1: # save MSA
				# save first seq as reprseq (fasta id)
				fprofile.write('>' + reprid + '\n')
				fprofile.write(str(clust_msa_seqs[nodes[0]]) + '\n')
				for node in nodes[1:]: # save all other seqs from cluster
					fprofile.write('>' + node + '\n')
					fprofile.write(str(clust_msa_seqs[node]) + '\n')
			else:
				# get singleton seq and save
				fprofile.write('>' + reprid + '\n')
				fprofile.write(str(cds_all[nodes[0]]) + '\n')
