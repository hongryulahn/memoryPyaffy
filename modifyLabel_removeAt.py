#!/usr/bin/env python
import sys
import argparse

parser=argparse.ArgumentParser(
        usage='''\
%(prog)s [options] exp sample2label out.txt

example: %(prog)s exp sample2label out.txt
''')

parser.add_argument('exp', help='expression file')
parser.add_argument('sample2label', help='sample ID file')
parser.add_argument('outfile', help='outfile')
args=parser.parse_args()


lst_nid = []
dic_nid2oid={}
IF=open(args.sample2label,'r')
for i, line in enumerate(IF):
	s=line.strip().split('\t')
	oid, nid = s[0].split('.')[0], s[1]
	lst_nid.append(nid)
	dic_nid2oid[nid]=oid

IF=open(args.exp,'r')
line = IF.readline()
s=line.strip().split('\t')
dic_oid2idx={}
for i, ss in enumerate(s[1:]):
	oid = ss.split('.')[0]
	if oid.startswith('GSM'):
		oid = oid.split('_')[0]
	dic_oid2idx[oid]=i

OF=open(args.outfile,'w')
OF.write('\t'.join([s[0]]+lst_nid)+'\n')

for line in IF:
	s=line.strip().split('\t')
	gene, lst_exp = s[0], s[1:]
	if gene.startswith('AFFX'):
		continue
	gene=gene.rstrip('_at')
	lst_newexp=[lst_exp[dic_oid2idx[dic_nid2oid[nid]]] for nid in lst_nid]
	OF.write('\t'.join([gene]+lst_newexp)+'\n')
