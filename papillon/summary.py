#!/usr/bin/env python
'''
  summarise papillon output
'''

import argparse
import collections
import csv
import logging
import sys

import numpy as np

def summary(in_fh, out_fh, metric, summary_fn):
  logging.info('starting...')

  coverage = collections.defaultdict(float) # average coverage
  size = collections.defaultdict(int) # average coverage
  for row in csv.DictReader(in_fh, delimiter='\t'):
    key = (row['Sample'], row['Gene'])
    size[key] += int(row['End']) - int(row['Start'])
    coverage[key] += (int(row['End']) - int(row['Start'])) * float(row[metric])

  # now write summary
  out = csv.DictWriter(out_fh, delimiter='\t', fieldnames=['Sample', 'Gene', metric])
  out.writeheader()
  genes = collections.defaultdict(list)
  for key in coverage:
    out.writerow({'Sample': key[0], 'Gene': key[1], metric: '{:.3f}'.format(coverage[key] / size[key])})
    genes[key[1]].append(coverage[key] / size[key])

  if summary_fn is not None:
    out = csv.DictWriter(open(summary_fn, 'w'), delimiter='\t', fieldnames=['Gene', 'n', metric, 'Min', 'Max', 'SD'])
    out.writeheader()
    for key in genes:
      out.writerow({'Gene': key, 'n': len(genes[key]), metric: np.mean(genes[key]), 'Min': min(genes[key]), 'Max': max(genes[key]), 'SD': np.std(genes[key])})

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Summary of papillon output by sample')
  parser.add_argument('--metric', required=False, default='Mean', help='which field to summarise')
  parser.add_argument('--summary', required=False, help='write overall summary')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  summary(sys.stdin, sys.stdout, args.metric, args.summary)
