#!/usr/bin/env python
'''
  generates z-values from papillon output
'''

import argparse
import collections
import csv
import logging
import sys

import numpy
import scipy.stats

def main(threshold):
  logging.info('reading from stdin...')

  inp = csv.DictReader(sys.stdin, delimiter='\t')
  groups = collections.defaultdict(list)
  rows = []
  for row in inp:
    # Sample  Chr     Start   End     Gene    Mean    Median  Min     Max     Pct
    key = (row['Chr'], row['Start'], row['End'])
    groups[key].append(float(row['Pct']))
    rows.append(row)

  logging.info('writing to stdout...')
  f = inp.fieldnames.copy()
  f.remove('Sample')
  f.remove('Gene')
  out = csv.DictWriter(sys.stdout, ['Sample', 'Gene'] + ['PctPValue', 'EventType', 'PctMean', 'PctSD', 'PctZ', 'PropMean'] + f, delimiter='\t')
  out.writeheader()
  for row in rows:
    key = (row['Chr'], row['Start'], row['End'])
    mean = sum(groups[key]) / len(groups[key])
    sd = numpy.std(groups[key])
    z = (float(row['Pct']) - mean) / sd
    p = scipy.stats.norm.sf(abs(z))*2
    if p < threshold:
      row['PctMean'] = mean
      row['PctSD'] = sd
      row['PctZ'] = z
      row['PctPValue'] = p
      row['PropMean'] = float(row['Pct']) / mean
      row['EventType'] = 'Deletion' if z < 0 else 'Duplication'
      out.writerow(row)

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Outlier assessment')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--significance_threshold', type=float, default=0.05, help='show with p-value below threshold')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.significance_threshold)
