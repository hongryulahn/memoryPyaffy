#!/usr/bin/env python
import sys
import argparse
import pyaffy
import datetime
import glob
import os
import logging
import collections
import itertools
import numpy as np
import time
from multiprocessing import Pool, Lock

logger = logging.getLogger(__name__)

def load_celfile(x):
	sample, cel_file, pm_idx, bg_correct, quantile_normalize = x
	isError=False
	try:
		y = pyaffy.celparser.parse_cel(cel_file)
		if sum(np.isnan(y)) != 0:
			isError=True
		else:
			mu = np.mean(y)
			if np.isinf(mu) or np.isnan(mu):
				isError=True
			std = np.std(y)
			if np.isinf(std) or np.isnan(std):
				isError=True
		if pm_idx is not None:
			y=y[pm_idx]
	except:
		isError=True
	if isError:
		return [sample, cel_file, isError, None, None, None]

	if bg_correct:
		Y=np.empty((y.size,1),dtype=np.float32)
		Y[:,0]=y
		Y=pyaffy.background.rma_bg_correct(Y)
		y=Y[:,0]

	if quantile_normalize:
		argsort_y=np.argsort(y)
		reverse_argsort_y=np.argsort(argsort_y)
	else:
		argsort_y=None

	return [sample, cel_file, isError, y, argsort_y, reverse_argsort_y]

def myBinaryColomnReader(IF, matrix_size, iter_nColumn, dtype, bufferSize=8*(2**30)):
	file_row_size, file_column_size = matrix_size
	nbytes=dtype(0).nbytes

	# define buffer_column_size and Y
	buffer_column_size = max(bufferSize/nbytes/file_row_size, 1)

	# initialize Y by reading from IF
	if file_column_size <= buffer_column_size and False:
		buffer_column_size=file_column_size
		Y=np.fromfile(IF, dtype=dtype).reshape(matrix_size)
	else:
		Y = np.empty((file_row_size, buffer_column_size), dtype = dtype)
		reading_column_size = min(file_column_size, buffer_column_size)
		IF.seek(0,0)
		for i in range(file_row_size):
			Y[i,0:reading_column_size] = np.fromstring(IF.read(reading_column_size*nbytes),dtype=dtype)
			IF.seek((file_column_size-reading_column_size)*nbytes,1)

	buffer_cur = 0
	file_cur = 0
	while True:
		nColumn = iter_nColumn.next()
		if file_cur + nColumn > file_column_size:
			print 'error: exceed entire column'
			sys.exit()
		if buffer_cur + nColumn <= buffer_column_size:
			buffer_cur += nColumn
			file_cur += nColumn
			yield Y[:,buffer_cur-nColumn:buffer_cur]
		else:
			# resize Y when nColumn exceeds buffer_column_size
			if nColumn > buffer_column_size:
				Y=np.empty((file_row_size, nColumn), dtype = dtype)
				buffer_column_size=nColumn

			if file_cur + buffer_column_size <= file_column_size:
				reading_column_size = buffer_column_size
			else:
				reading_column_size = file_column_size - file_cur
			IF.seek(file_cur*nbytes,0)
			for i in range(file_row_size):
				Y[i,0:reading_column_size] = np.fromstring(IF.read(reading_column_size*nbytes),dtype=dtype)
				IF.seek((file_column_size-reading_column_size)*nbytes,1)
			buffer_cur = 0
	
def myRMA(
	cdf_file,
	lst_sample,
	lst_celFile,
	error_dir,
	output_file,
	pm_probes_only = True,
	bg_correct = True,
	quantile_normalize = True,
	log2_normalize = True,
	medianpolish = True,
	MEMSIZE = 30,
	dtype=np.float32,
	nProcess = 1
	):
	"""Perform RMA on a set of samples.

	Parameters
	----------
	cdf_file: str or unicode
		The path of the Brainarray CDF file to use.
		Note: Brainarray CDF files can be downloaded from
			http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/genomic_curated_CDF.asp
	lst_sample: sample name
	lst_celFile: CEL file (must be gzip'ed) corresponded to lst_sample
	error_dir:  directory where the error file moves to
	output_file: output file
	pm_probes_only: bool, optional
		Whether or not to only use PM (perfect match) probes and ignore all MM
		(mismatch) probes. [True]
	bg_correct: bool, optional
		Whether or not to apply background correction. [True]
	quantile_normalize: bool, optional
		Whether or not to apply quantile normalization. [True]
	log2_normalize: bool, optional
		Whether or not to apply log2 normalization. [True]
	medianpolish: bool, optional
		Whether or not to apply medianpolish. [True]
	MEMSIZE: the limit of memory size to be used (Gigabytes) [30]
	dtype:  datatype that numpy array uses [numpy.float32]
	nProcess: int, optional
		The number of process for multiprocessing

	Returns
	-------
	

	Examples
	--------
	>>> from collections import OrderedDict
	>>> import pyaffy
	>>> cdf_file = '/path/to/brainarray/cdf/HGU133Plus2_Hs_ENTREZG.cdf'
	>>> sample_cel_files = OrderedDict([
			['Sample 1', '/path/to/sample_1.CEL.gz'],
			['Sample 2', '/path/to/sample_2.CEL.gz'],
		])
	>>> genes, samples, X = pyaffy.rma(cdf_file, sample_cel_files)
	"""

	### checks
	assert (isinstance(cdf_file, (str, unicode)) and os.path.isfile(cdf_file)) or (isinstance(cdf_file, tuple) and len(cdf_file)==4)

	assert isinstance(lst_sample, list)
	for sample in lst_sample:
		assert isinstance(sample, (str, unicode))

	assert isinstance(lst_celFile, list)
	for celFile in lst_celFile:
		assert isinstance(celFile, (str, unicode))
		assert os.path.isfile(celFile), \
				'CEL file "%s" does not exist!' %(cel_file)

	assert isinstance(error_dir, (str, unicode))
	assert isinstance(output_file, (str, unicode))
	assert isinstance(nProcess, int)

	# MEMSIZE transition into Gigabytes (1G < MEMSIZE < 1000G)
	MEMSIZE = min(max(MEMSIZE, 1),1000)*2**30

	t00 = time.time()

	### read CDF data
	logger.info('BEGINNING: Parsing CDF file.')
	t0 = time.time()
	# parse the CDF file
	if isinstance(cdf_file, (str, unicode)):
		probe_type = 'pm'
		if not pm_probes_only:
			probe_type = 'all'
		name, num_rows, num_cols, pm_probesets = \
				pyaffy.cdfparser.parse_cdf(cdf_file,probe_type = probe_type)
	elif isinstance(cdf_file, tuple) and len(cdf_file) == 4:
		name, num_rows, num_cols, pm_probesets = cdf_file
	pm_probesets = collections.OrderedDict(sorted(pm_probesets.items()))

	# concatenate indices of all PM probes into one long vector
	pm_sel = np.concatenate(pm_probesets.values())

	t1 = time.time()
	logger.info('CDF file parsing time: %.2f s', t1 - t0)
	logger.info('CDF array design name: %s', name)
	logger.info('CDF rows / columns: %d x %d', num_rows, num_cols)
	logger.info('Using only use PM (perfect match) probes: %s'%(pm_probes_only))
	logger.info('The number of PM probes: %s'%(len(pm_sel)))
	logger.info('ENDING: Parsing CDF file.')

	### read CEL data
	logger.info('BEGINNING: Parsing %d CEL files...'%(len(lst_celFile)))
	t0 = time.time()

	if nProcess > 1:
		pool = Pool(processes = nProcess)
		lst_result = pool.imap(load_celfile, itertools.izip(iter(lst_sample),iter(lst_celFile),iter([pm_sel]*len(lst_sample)),iter([True]*len(lst_sample)),iter([True]*len(lst_sample))), chunksize=10)
		pool.close()
	else:
		lst_result = itertools.imap(load_celfile, itertools.izip((lst_sample),iter(lst_celFile),iter([pm_sel]*len(lst_sample)),iter([True]*len(lst_sample)),iter([True]*len(lst_sample))))

	nProbe = pm_sel.size
	probematrix_file = output_file+'.probematrix.binary'
	OF=open(output_file+'.probematrix.binary', 'wb')
	if quantile_normalize:
		if np.uint32(nProbe) == nProbe:
			probematrix_dtype = np.uint32
		else:
			probematrix_dtype = np.uint64
		quantile_sum= np.zeros(nProbe, dtype = np.float32)
	else:
		probematrix_dtype = dtype
	lst_correct_sample=[]
	nCorrect,nError=0,0
	for j,result in enumerate(lst_result):
		sample, celFile, isError, y, argsort_y, reverse_argsort_y = result
		if isError:
			nError+=1
			if not os.path.isdir(error_dir):
				os.system('mkdir '+error_dir)
			os.system('mv -f '+celFile+' '+error_dir)
		else:
			lst_correct_sample.append(str(sample))
			if quantile_normalize:
				quantile_sum = np.add(quantile_sum, y[argsort_y])
				probematrix_dtype(reverse_argsort_y).tofile(OF)
			else:
				if log2_normalize:
					y = np.log2(y)
				probematrix_dtype(y).tofile(OF)
			del result[:]
			nCorrect+=1
		logger.info('Parsed %d/%d CEL file: filename=%s, isError=%s'%(j+1,len(lst_celFile),celFile,isError))
	OF.close()

	if quantile_normalize:
		quantile_mean = quantile_sum / np.float(nCorrect)
		if log2_normalize:
			quantile_mean = np.log2(quantile_mean)
	t1 = time.time()
	logger.info('CEL file parsing time: %.2f s', t1 - t0)
	logger.info('The number of CEL files: corrected=%d, error=%d, total tested=%d'%(nCorrect, nError, nCorrect+nError))
	logger.info('Background correction: %s'%(bg_correct))
	logger.info('Quantile calculation: %s'%(quantile_normalize))
	logger.info('The shape of temporary probe matrix: %s * %s'%(nCorrect, len(pm_sel)))
	logger.info('Temporary probe matrix file: %s'%(probematrix_file))
	logger.info('ENDING: Parsing %d CEL files...'%(len(lst_celFile)))

	### probeset summarization (with or without median polish)
	logger.info('BEGINNING: Summarizing Probes into Genes or ProbeSets...')

	#return map(str,pm_sel), samples, Y
	t0 = time.time()
	genes=pm_probesets.keys()

	OF=open(output_file, 'w')
	OF.write('IDs\t'+'\t'.join(lst_correct_sample)+'\n')

	num_converged = 0
	IF=open(probematrix_file, 'rb')
	for gidx, X_sub in enumerate(myBinaryColomnReader(IF, matrix_size=(len(lst_correct_sample),nProbe), iter_nColumn=(probes.size for probes in pm_probesets.itervalues()), dtype=probematrix_dtype, bufferSize=MEMSIZE)):
		if quantile_normalize:
			Y_sub = np.empty(X_sub.shape, dtype=np.float32)
			for i in range(X_sub.shape[0]):
				Y_sub[i,:] = quantile_mean[X_sub[i,:]]
		else:
			Y_sub = X_sub
		Y_sub = np.transpose(Y_sub)
		if medianpolish:
			_, row_eff, col_eff, global_eff, converged, num_iter = pyaffy.medpolish.medpolish(Y_sub, copy = False)
			Y = col_eff + global_eff
			if converged:
				num_converged += 1
		else:
			# simply use median across probes
			Y = np.median(Y_sub, axis = 0)
		#if sorted_output:
		OF.write(genes[gidx]+'\t'+'\t'.join(map(str,Y))+'\n')
	OF.close()

	t1 = time.time()
	logger.info('Probeset summarization time: %.2f s.', t1 - t0)
	logger.info('Using medianpolish: %s'%(medianpolish))
	if medianpolish:
		logger.debug('Median polish Converged: %d / %d (%.1f%%)',
				num_converged, nProbe, 100 * (num_converged / float(nProbe)))
	logger.info('ENDING: Summarizing Probes into Genes or ProbeSets...')


	### report total time
	t11 = time.time()
	logger.info('Total RMA time: %.1f s.', t11 - t00)
	logger.info('Excution date: %s', str(datetime.date.today()))
	logger.info('OPTION: Using only use PM (perfect match) probes: %s'%(pm_probes_only))
	logger.info('OPTION: Background correction: %s'%(bg_correct))
	logger.info('OPTION: Quantile normalization: %s'%(quantile_normalize))
	logger.info('OPTION: Log2 scale normalization: %s'%(log2_normalize))
	logger.info('OPTION: Using medianpolish summarization: %s'%(medianpolish))
	logger.info('The number of IDs (Genes or Probe sets): %d'%(len(genes)))
	logger.info('The number of CEL files: corrected=%d, error=%d, total tested=%d'%(nCorrect, nError, nCorrect+nError))
	logger.info('The shape of output matrix: %s * %s'%(len(genes), nCorrect))
	logger.info('output matrix file: %s'%(output_file))
	logger.info('ENDING: Total RMA')


if __name__ == '__main__':
	parser=argparse.ArgumentParser(
		usage='''\
%(prog)s [options] rice.cdf raw_data out.txt

example: %(prog)s rice.cdf raw_data out.txt -p 10
''')
	parser.add_argument('cdf', help='cdf file')
	parser.add_argument('dir', help='raw data directory')
	parser.add_argument('outfile', help='outfile')
	parser.add_argument('-log', required=False, metavar='str', default='None', help='log file')
	parser.add_argument('-MEMSIZE', type=int, required=False, metavar='N', default='8', help='limit of memory usage')
	parser.add_argument('-p', type=int, required=False, metavar='N', default='1', help='total process')
	args=parser.parse_args()

	if args.log == 'None':
		logfile=args.outfile+'.log'
	else:
		logfile=args.log
	logging.basicConfig(filename=logfile, filemode='w', level=logging.INFO)

	lst_celFile = sorted(list(glob.glob(args.dir+'/*CEL*'))+list(glob.glob(args.dir+'/*cel*')))
	lst_sample = [os.path.basename(celFile) for celFile in lst_celFile]
	error_dir = args.dir+'/error/'
	output_file = args.outfile

	myRMA(
		cdf_file=args.cdf,
		lst_sample=lst_sample,
		lst_celFile=lst_celFile,
		error_dir=error_dir,
		output_file=output_file,
		pm_probes_only = True,
		bg_correct = True,
		quantile_normalize = True,
		log2_normalize = True,
		medianpolish = True,
		MEMSIZE = args.MEMSIZE,
		dtype=np.float32,
		nProcess = args.p 
	)
