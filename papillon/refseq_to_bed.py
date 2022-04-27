#!/usr/bin/env python
'''
  list of genes to exons from refseq
  e.g.
  python refseq_to_bed.py --genes MSH2 RNF43 --refseq ~/crc/data/public_datasets/refseq.ucsc.hg19.180829.gz
  #bin    name    chrom   strand  txStart txEnd   cdsStart        cdsEnd  exonCount       exonStarts      exonEnds        score   name2   cdsStartStat    cdsEndStat      exonFrames
  0       NM_032291       chr1    +       66999638        67216822        67000041        67208778        25      66999638,67091529,67098752,67101626,67105459,67108492,67109226,67126195,67133212,67136677,67137626,67138963,67142686,67145360,67147551,67154830,67155872,67161116,67184976,67194946,67199430,67205017,67206340,67206954,67208755,     67000051,67091593,67098777,67101698,67105516,67108547,67109402,67126207,67133224,67136702,67137678,67139049,67142779,67145435,67148052,67154958,67155999,67161176,67185088,67195102,67199563,67205220,67206405,67207119,67216822,     0       SGIP1   cmpl    cmpl    0,1,2,0,0,0,1,0,0,0,1,2,1,1,1,1,0,1,1,2,2,0,2,1,1,
'''

import argparse
import collections
import csv
import gzip
import logging
import sys

def main(genes, transcripts, refseq, coding):
  logging.info('starting...')
  genes = set(genes)
  if transcripts is not None:
    transcripts = set(transcripts)
  mins = collections.defaultdict(int)
  mins.default_factory = lambda: 3e9
  maxs = collections.defaultdict(int)
  for row in csv.DictReader(gzip.open(refseq, 'rt'), delimiter='\t'):
    if (genes is not None or row['name2'] in genes) and (transcripts is None or row['name'] in transcripts):
      for x, y in zip(row['exonStarts'].split(','), row['exonEnds'].split(',')):
        if x != '' and y != '':
          if coding:
            min_coding = int(row['cdsStart'])
            max_coding = int(row['cdsEnd'])
            x = str(max([int(x), min_coding]))
            y = str(min([int(y), max_coding]))
          sys.stdout.write('{}\t{}\t{}\t{}\n'.format(row['chrom'], x, y, row['name2']))

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Generate region for specified genes')
  parser.add_argument('--genes', required=False, nargs='+', help='genes to include')
  parser.add_argument('--transcripts', required=False, nargs='+', help='transcripts to include')
  parser.add_argument('--refseq', required=True, help='refseq gz')
  parser.add_argument('--coding', action='store_true', help='coding regions only')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.genes, args.transcripts, args.refseq, args.coding)
