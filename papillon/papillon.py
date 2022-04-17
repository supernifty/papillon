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

def statistic(pileups, stat, max_value=None, mapped=None):
  if stat == 'mean':
    mean = sum(pileups) / len(pileups)
    result = mean
  elif stat == 'min':
    result = min(pileups)
  elif stat == 'max':
    result = max(pileups)
  elif stat == 'median':
    if len(pileups) % 2 == 0:
      median = (pileups[int(len(pileups) / 2)] + pileups[int(len(pileups) / 2) - 1]) / 2
    else:
      median = pileups[int(len(pileups) / 2)]
    result = median
  elif stat == 'percent':
    result = sum(pileups) * 100 / mapped # TODO we don't consider read length
  elif stat == 'sd':
    result = np.std(pileups, ddof=1)

  if max_value is not None:
    return min(result, max_value)
  else:
    return result

def main(bams, bed, genes_list, plot, capture, stat, exon_plots, padding, max_coverage, min_mapq, min_coverage, base_level, sample_level, raw, sample_name_end):
  logging.info('starting...')

  if genes_list is not None:
    genes = set(genes_list)
    genes_list = sorted(genes_list)
    logging.info('limiting to %i genes', len(genes))
  else:
    genes = None # no limit

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
    if len(fields) < 3:
      continue
    bed_start = int(fields[1]) - padding
    bed_finish = int(fields[2]) + padding
    if capture is not None and fields[0] in captures:
      intersections = captures[fields[0]][bed_start:bed_finish]
      total = 0
      for interval in intersections:
        total += min(interval.end, bed_finish) - max(interval.begin, bed_start)
    else:
      overlap = ''
    if len(fields) < 4:
      fields.append('none')
    regions[fields[3]].add((fields[0], bed_start, bed_finish, overlap))
  logging.info('processing %s: done %i lines found %i genes', bed, lines, len(regions))

  # select all genes
  if genes is None:
    genes = set(regions.keys())
    genes_list = sorted(list(genes))

  if plot is not None:
    gene_result = {}
    gene_result_annot = {}
    combined_result = np.zeros((len(bams), len(genes)), dtype=float)
    for gene in regions:
      gene_result[gene] = np.zeros((len(bams), len(regions[gene])), dtype=float)
      gene_result_annot[gene] = np.empty((len(bams), len(regions[gene])), dtype='S5')

  if base_level is not None:
    bases = collections.defaultdict(list)
    base_fh = open(base_level, 'w')
    base_fh.write('Chr\tPos\tn\tMean\tMedian\tMin\tMax\tSD\tpct2_5\tpct97_5\n')

  if raw is not None:
    raw_fh = open(raw, 'w')
    raw_fh.write('Sample\tChr\tPos\tDepth\n')

  sys.stdout.write('Sample\tChr\tStart\tEnd\tGene\tMean\tMedian\tMin\tMax\tPct\n')
  
  xticklabels = collections.defaultdict(list) # genes to list of exons
  yticklabels = []

  # process each bam
  bam_idx = 0
  for bam_idx, bam in enumerate(bams):
    logging.info('processing file %i of %i: %s...', bam_idx + 1, len(bams), bam)
    samfile = pysam.AlignmentFile(bam, "rb" )
    sample_name = bam.split('/')[-1].split(sample_name_end)[0]
    yticklabels.append(sample_name)
    gene_count = 0
    gene_pileups = collections.defaultdict(list)
    for gene_count, gene in enumerate(regions):
      if gene in genes:
        for region_idx, region in enumerate(sorted(regions[gene], key=lambda x: x[1])): # each region in the gene
          if bam_idx == 0:
            if plot is not None:
              xticklabels[gene].append('{}\n{}bp{}'.format(region[1], region[2]-region[1], region[3]))
          # note: pysam by default filters duplicates
          pileups = [x.n for x in samfile.pileup(region[0], region[1], region[2], min_mapping_quality=min_mapq) if region[1] <= x.pos < region[2]]
          logging.debug(pileups)
          if len(pileups) < (region[2] - region[1]):
            pileups += [0] * (region[2] - region[1])
          gene_pileups[gene] += pileups
          if base_level is not None or raw is not None:
            for i, r in enumerate(range(region[1], region[2])):
              if base_level is not None:
                bases[(region[0], r)].append(pileups[i])
              if raw is not None:
                raw_fh.write('{}\t{}\t{}\t{}\n'.format(bam, region[0], r, pileups[i]))
          sorted_pileups = sorted(pileups)
          median = statistic(sorted_pileups, 'median')
          mean = statistic(sorted_pileups, 'mean')
          percent_of_mapped = statistic(pileups, 'percent', None, samfile.mapped)
          sys.stdout.write('{}\t{}\t{}\t{}\t{}\t{:.1f}\t{:.1f}\t{}\t{}\t{:.6f}\n'.format(bam, region[0], region[1], region[2], gene, mean, median, min(pileups), max(pileups), percent_of_mapped))
          if plot is not None:
            gene_result[gene][bam_idx, region_idx] = statistic(sorted_pileups, stat, max_coverage, samfile.mapped)
            if max_coverage is not None and gene_result[gene][bam_idx, region_idx] >= max_coverage:
              gene_result_annot[gene][bam_idx, region_idx] = '{}+'.format(max_coverage)
            elif min_coverage is not None and gene_result[gene][bam_idx, region_idx] >= min_coverage:
              gene_result_annot[gene][bam_idx, region_idx] = ''
            else:
              gene_result_annot[gene][bam_idx, region_idx] = '{:.0f}'.format(gene_result[gene][bam_idx, region_idx])
      if gene_count % 1000 == 0:
        logging.debug('%i genes processed', gene_count + 1)
    logging.info('%i genes processed', gene_count + 1)

    # now deal with gene pileup
    if plot is not None:
      for gene in gene_pileups:
        gene_pileup = sorted(gene_pileups[gene])
        combined_result[bam_idx, genes_list.index(gene)] = statistic(gene_pileup, stat, max_coverage, samfile.mapped)
        logging.debug('updated sample %i %s %s', bam_idx, sample_name, gene)

  if base_level is not None:
    for pos in bases:
      base_fh.write('{}\t{}\t{}\t{:.1f}\t{}\t{}\t{}\t{:.1f}\t{:.1f}\t{:.1f}\n'.format(pos[0], pos[1], len(bases[pos]), statistic(bases[pos], 'mean'), statistic(bases[pos], 'median'), statistic(bases[pos], 'min'), statistic(bases[pos], 'max'), statistic(bases[pos], 'sd'), np.percentile(bases[pos], 2.5), np.percentile(bases[pos], 97.5)))

  if sample_level is not None:
    logging.info('writing to %s...', sample_level)
    base_len = len(bases)
    with open(sample_level, 'w') as sl_fh:
      sl_fh.write('Sample\tBases\tMean\n')
      for i in range(0, bam_idx + 1):
        total = sum([bases[x][i] for x in bases])
        sl_fh.write('{}\t{}\t{:.2f}\n'.format(bams[i], base_len, total / base_len))

  if plot is not None:
    logging.info('plotting...')

    if exon_plots:
      # each gene individually
      for gene in regions:
        if genes is None or gene in genes:
          target_image = '{}.{}.png'.format(plot, gene)
          logging.info('plotting %s with %i x %i...', target_image, gene_result[gene].shape[1], gene_result[gene].shape[0])
          fig, ax = plt.subplots(figsize=(max(FIGSIZE_FACTOR_X_MIN + FIGSIZE_FACTOR_X * 6, FIGSIZE_FACTOR_X_MIN + FIGSIZE_FACTOR_X * gene_result[gene].shape[1] / 4), max(FIGSIZE_FACTOR_Y * 8, FIGSIZE_FACTOR_Y * gene_result[gene].shape[0] / 10)))
          extra = ''
          if padding is not None:
            extra += 'Padding {}bp. '.format(padding)
          if max_coverage is not None:
            extra += 'Max coverage {}. '.format(max_coverage)
          ax.set_title('Coverage plot for {} with {} region coverage. {}'.format(gene, stat, extra))
          heatmap = sns.heatmap(gene_result[gene], xticklabels=xticklabels[gene], yticklabels=yticklabels, annot=np.char.decode(gene_result_annot[gene]), ax=ax, cmap="RdYlGn", fmt='', vmin=0)
          ax.set_xlabel('Regions') # TODO doesn't work
          ax.set_ylabel('Samples') # TODO doesn't work
          fig = heatmap.get_figure()
          fig.savefig(target_image)

    # all genes
    target_image = '{}.png'.format(plot)
    logging.info('plotting %s with %i x %i: %s %s...', target_image, combined_result.shape[1], combined_result.shape[0], genes_list, yticklabels)
    fig, ax = plt.subplots(figsize=(max(FIGSIZE_FACTOR_X_MIN + FIGSIZE_FACTOR_X * 6, FIGSIZE_FACTOR_X_MIN + FIGSIZE_FACTOR_X * combined_result.shape[1] / 4), max(FIGSIZE_FACTOR_Y * 8, FIGSIZE_FACTOR_Y * combined_result.shape[0] / 10)))
    ax.set_title('Coverage plot with {} region coverage'.format(stat))
    heatmap = sns.heatmap(combined_result, xticklabels=genes_list, yticklabels=yticklabels, annot=True, ax=ax, fmt='.0f', cmap="RdYlGn", vmin=0)
    ax.set_xlabel('Regions') # TODO doesn't work
    ax.set_ylabel('Samples') # TODO doesn't work
    fig = heatmap.get_figure()
    fig.savefig(target_image)
 
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Measure coverage across exons and genes')
  parser.add_argument('--bams', required=True, nargs='+', help='bams to analyse')
  parser.add_argument('--genes', required=False, nargs='*', help='genes to filter on')
  parser.add_argument('--bed', required=True, help='regions of interest')
  parser.add_argument('--stat', required=False, default='mean', choices=('mean', 'min', 'max', 'median', 'percent'), help='value to report on plot')
  parser.add_argument('--capture', required=False, help='capture to compare to')
  parser.add_argument('--plot', required=False, help='graph file prefix e.g. heatmap will generate prefix.GENE.png')
  parser.add_argument('--padding', required=False, default=0, type=int, help='padding to apply to bed')
  parser.add_argument('--min_mapq', required=False, default=0, type=int, help='minimum allowable mapping quality for read')
  parser.add_argument('--max_coverage', required=False, default=None, type=int, help='do not report coverage above this on plots')
  parser.add_argument('--min_coverage', required=False, default=None, type=int, help='leave empty value if above this')
  parser.add_argument('--exon_plots', action='store_true', help='include exon plots')
  parser.add_argument('--base_level', required=False, help='filename for base level coverage')
  parser.add_argument('--sample_level', required=False, help='filename for sample level coverage')
  parser.add_argument('--sample_name_end', required=False, default='.', help='string marking end of sample name')
  parser.add_argument('--raw', required=False, help='filename for individual base/sample coverage')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.bams, args.bed, args.genes, args.plot, args.capture, args.stat, args.exon_plots, args.padding, args.max_coverage, args.min_mapq, args.min_coverage, args.base_level, args.sample_level, args.raw, args.sample_name_end)
