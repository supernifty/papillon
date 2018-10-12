#!/usr/bin/env python
'''
  measure coverage across regions and bams
'''

import argparse
import collections
import logging
import sys

import numpy as np
import pysam
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

FIGSIZE_FACTOR_X=3.0
FIGSIZE_FACTOR_Y=1.0

def main(bams, bed, genes_list, plot):
  logging.info('starting...')

  if genes_list is not None:
    genes = set(genes_list)
    logging.info('limiting to %i genes', len(genes))
  else:
    genes = None

  regions = collections.defaultdict(set)
  logging.info('processing %s...', bed)
  lines = 0
  for lines, line in enumerate(open(bed, 'r')):
    if line.startswith('#'):
      continue
    fields = line.strip('\n').split('\t')
    regions[fields[3]].add((fields[0], int(fields[1]), int(fields[2])))
  logging.info('processing %s: done %i lines found %i genes', bed, lines, len(regions))

  if plot is not None:
    gene_result = {}
    for gene in regions:
      gene_result[gene] = np.zeros((len(bams), len(regions[gene])), dtype=float)

  sys.stdout.write('Sample\tChr\tStart\tEnd\tGene\tMean\tMedian\tMin\tMax\n')
  
  xticklabels = collections.defaultdict(list)
  yticklabels = []
  for bam_idx, bam in enumerate(bams):
    logging.info('processing %s', bam)
    samfile = pysam.AlignmentFile(bam, "rb" )
    sample_name = bam.split('/')[-1].split('.')[0]
    yticklabels.append(sample_name)
    for gene_count, gene in enumerate(regions):
      if genes is None or gene in genes:
        for region_idx, region in enumerate(sorted(regions[gene], key=lambda x: x[1])): # each region
          if bam_idx == 0:
            xticklabels[gene].append(region[1])
          pileups = [x.n for x in samfile.pileup(region[0], region[1], region[2]) if region[1] <= x.pos < region[2]]
          logging.debug(pileups)
          if len(pileups) < (region[2] - region[1]):
            pileups += [0] * (region[2] - region[1])
          sorted_pileups = sorted(pileups)
          if len(pileups) % 2 == 0:
            median = (sorted_pileups[int(len(pileups) / 2)] + sorted_pileups[int(len(pileups) / 2) - 1]) / 2
          else:
            median = sorted_pileups[int(len(pileups) / 2)]
          sys.stdout.write('{}\t{}\t{}\t{}\t{}\t{:.1f}\t{:.1f}\t{}\t{}\n'.format(bam, region[0], region[1], region[2], gene, sum(pileups) / len(pileups), median, min(pileups), max(pileups)))
          if plot is not None:
            gene_result[gene][bam_idx, region_idx] = int(median)
      if gene_count % 100 == 0:
        logging.info('%i genes processed', gene_count)

  if plot is not None:
    logging.info('plotting...')
    for gene in regions:
      if genes is None or gene in genes:
        target_image = '{}.{}.png'.format(plot, gene)
        logging.info('plotting %s with %i x %i...', target_image, gene_result[gene].shape[1], gene_result[gene].shape[0])
        fig, ax = plt.subplots(figsize=(max(FIGSIZE_FACTOR_X * 6, FIGSIZE_FACTOR_X * gene_result[gene].shape[1] / 10), max(FIGSIZE_FACTOR_Y * 8, FIGSIZE_FACTOR_Y * gene_result[gene].shape[0] / 10)))
        heatmap = sns.heatmap(gene_result[gene], xticklabels=xticklabels[gene], yticklabels=yticklabels, annot=True, ax=ax, fmt='.0f', cmap="coolwarm_r", vmin=0)
        fig = heatmap.get_figure()
        fig.savefig(target_image)

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--bams', required=True, nargs='+', help='bams to analyse')
  parser.add_argument('--genes', required=False, nargs='*', help='genes to filter on')
  parser.add_argument('--bed', required=True, help='regions of interest')
  parser.add_argument('--plot', required=False, help='graph file prefix e.g. heatmap will generate heatmap.GENE.png')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.bams, args.bed, args.genes, args.plot)
