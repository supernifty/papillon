#!/usr/bin/env python
'''
  measure coverage across regions and bams
'''

import argparse
import collections
import logging
import sys

import intervaltree
import numpy as np
import pysam
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

FIGSIZE_FACTOR_X=3.0
FIGSIZE_FACTOR_X_MIN=0

FIGSIZE_FACTOR_Y=1.5

def statistic(pileups, stat):
  if stat == 'mean':
    mean = sum(pileups) / len(pileups)
    return int(mean)
  elif stat == 'min':
    return min(pileups)
  elif stat == 'max':
    return max(pileups)
  elif stat == 'median':
    if len(pileups) % 2 == 0:
      median = (pileups[int(len(pileups) / 2)] + pileups[int(len(pileups) / 2) - 1]) / 2
    else:
      median = pileups[int(len(pileups) / 2)]
    return int(median)

def main(bams, bed, genes_list, plot, capture, stat):
  logging.info('starting...')

  if genes_list is not None:
    genes = set(genes_list)
    logging.info('limiting to %i genes', len(genes))
  else:
    genes = None

  captures = {}
  if capture is not None:
    capture_line = 0
    logging.info('processing %s...', capture)
    for capture_line, line in enumerate(open(capture, 'r')):
      if line.startswith('#'):
        continue
      fields = line.strip('\n').split('\t')
      if len(fields) < 3:
        continue
      if fields[0] not in captures:
        captures[fields[0]] = intervaltree.IntervalTree()
      captures[fields[0]][int(fields[1]):int(fields[2])] = True
    logging.info('%i capture lines processed', capture_line)  

  regions = collections.defaultdict(set)
  logging.info('processing %s...', bed)
  lines = 0
  for lines, line in enumerate(open(bed, 'r')):
    if line.startswith('#'):
      continue
    fields = line.strip('\n').split('\t')
    if capture is not None and fields[0] in captures:
      intersections = captures[fields[0]][int(fields[1]):int(fields[2])]
      total = 0
      for interval in intersections:
        total += min(interval.end, int(fields[2])) - max(interval.begin, int(fields[1]))
      overlap = '\n{:.0f}%'.format(100. * total / (int(fields[2]) - int(fields[1])))
    else:
      overlap = ''
    regions[fields[3]].add((fields[0], int(fields[1]), int(fields[2]), overlap))
  logging.info('processing %s: done %i lines found %i genes', bed, lines, len(regions))

  if plot is not None:
    gene_result = {}
    combined_result = np.zeros((len(bams), len(genes)), dtype=float)
    for gene in regions:
      gene_result[gene] = np.zeros((len(bams), len(regions[gene])), dtype=float)

  sys.stdout.write('Sample\tChr\tStart\tEnd\tGene\tMean\tMedian\tMin\tMax\n')
  
  xticklabels = collections.defaultdict(list) # genes to list of exons
  yticklabels = []
  for bam_idx, bam in enumerate(bams):
    logging.info('processing file %i of %i: %s...', bam_idx + 1, len(bams), bam)
    samfile = pysam.AlignmentFile(bam, "rb" )
    sample_name = bam.split('/')[-1].split('.')[0]
    yticklabels.append(sample_name)
    gene_count = 0
    for gene_count, gene in enumerate(regions):
      if genes is None or gene in genes:
        gene_pileup = []
        for region_idx, region in enumerate(sorted(regions[gene], key=lambda x: x[1])): # each region in the gene
          if bam_idx == 0:
            xticklabels[gene].append('{}\n{}bp{}'.format(region[1], region[2]-region[1], region[3]))
          pileups = [x.n for x in samfile.pileup(region[0], region[1], region[2]) if region[1] <= x.pos < region[2]]
          logging.debug(pileups)
          if len(pileups) < (region[2] - region[1]):
            pileups += [0] * (region[2] - region[1])
          gene_pileup += pileups
          sorted_pileups = sorted(pileups)
          median = statistic(sorted_pileups, 'median')
          mean = statistic(sorted_pileups, 'mean')
          sys.stdout.write('{}\t{}\t{}\t{}\t{}\t{:.1f}\t{:.1f}\t{}\t{}\n'.format(bam, region[0], region[1], region[2], gene, mean, median, min(pileups), max(pileups)))
          if plot is not None:
            gene_result[gene][bam_idx, region_idx] = statistic(sorted_pileups, stat)
        # now deal with gene pileup
        gene_pileup = sorted(gene_pileup)
        combined_result[bam_idx, genes_list.index(gene)] = statistic(gene_pileup, stat)

      if gene_count % 1000 == 0:
        logging.info('%i genes processed', gene_count)
    logging.info('%i genes processed', gene_count)

  if plot is not None:
    logging.info('plotting...')

    # each gene individually
    for gene in regions:
      if genes is None or gene in genes:
        target_image = '{}.{}.png'.format(plot, gene)
        logging.info('plotting %s with %i x %i...', target_image, gene_result[gene].shape[1], gene_result[gene].shape[0])
        fig, ax = plt.subplots(figsize=(max(FIGSIZE_FACTOR_X_MIN + FIGSIZE_FACTOR_X * 6, FIGSIZE_FACTOR_X_MIN + FIGSIZE_FACTOR_X * gene_result[gene].shape[1] / 4), max(FIGSIZE_FACTOR_Y * 8, FIGSIZE_FACTOR_Y * gene_result[gene].shape[0] / 10)))
        ax.set_title('Coverage plot for {} with {} region coverage'.format(gene, stat))
        heatmap = sns.heatmap(gene_result[gene], xticklabels=xticklabels[gene], yticklabels=yticklabels, annot=True, ax=ax, fmt='.0f', cmap="Spectral", vmin=0)
        ax.set_xlabel('Regions') # TODO doesn't work
        ax.set_ylabel('Samples') # TODO doesn't work
        fig = heatmap.get_figure()
        fig.savefig(target_image)

    # all genes
    target_image = '{}.png'.format(plot)
    logging.info('plotting %s with %i x %i...', target_image, combined_result.shape[1], combined_result.shape[0])
    fig, ax = plt.subplots(figsize=(max(FIGSIZE_FACTOR_X_MIN + FIGSIZE_FACTOR_X * 6, FIGSIZE_FACTOR_X_MIN + FIGSIZE_FACTOR_X * combined_result.shape[1] / 4), max(FIGSIZE_FACTOR_Y * 8, FIGSIZE_FACTOR_Y * combined_result.shape[0] / 10)))
    ax.set_title('Coverage plot with {} region coverage'.format(stat))
    heatmap = sns.heatmap(combined_result, xticklabels=genes, yticklabels=yticklabels, annot=True, ax=ax, fmt='.0f', cmap="Spectral", vmin=0)
    ax.set_xlabel('Regions') # TODO doesn't work
    ax.set_ylabel('Samples') # TODO doesn't work
    fig = heatmap.get_figure()
    fig.savefig(target_image)
 
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--bams', required=True, nargs='+', help='bams to analyse')
  parser.add_argument('--genes', required=False, nargs='*', help='genes to filter on')
  parser.add_argument('--bed', required=True, help='regions of interest')
  parser.add_argument('--stat', required=False, default='mean', choices=('mean', 'min', 'max', 'median'), help='capture to compare to')
  parser.add_argument('--capture', required=False, help='capture to compare to')
  parser.add_argument('--plot', required=False, help='graph file prefix e.g. heatmap will generate prefix.GENE.png')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.bams, args.bed, args.genes, args.plot, args.capture, args.stat)
