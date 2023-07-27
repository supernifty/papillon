#!/usr/bin/env python
'''
  measure coverage across regions and bams
  - takes main papillon output as its input
'''

import argparse
import csv
import collections
import logging
import math
import sys

import intervaltree
import numpy as np
import pysam
import matplotlib as mpl
from matplotlib.ticker import FuncFormatter
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

MAX_LOSS=-5
CLASSIFY=(2, 3) # two out of three consecutive outside range
LOSS=-0.5
GAIN=0.5

LOSS_COLOR='red'
GAIN_COLOR='green'

#VALUE='Pct'
#VALUE='Mean'

def fmt_container(m):
  def x_fmt(x, pos):
    if x in m:
      return m[x]
    else:
      return ''
  return x_fmt

def main(target, title, highlight, dpi=300, width=12, height=8, xaxis=None, yaxis=None, chromosome=None, background=None, pct=False, normalise=False):
  logging.info('starting...')

  if pct:
    colname = 'Pct'
  else:
    colname = 'Mean'

  data = {'x': [], 'max': [], 'mean': [], 'min': []} # max, mean, min
  chromosomes = {}
  last_chrom = None
  chrom_offset = 0
  last_pos = 0
  fig, ax = plt.subplots()
  # not working
  #ax.xaxis.set_major_formatter(FuncFormatter(fmt_container(chromosomes)))
  #ax.xaxis.set_minor_formatter(FuncFormatter(fmt_container(chromosomes)))

  # need to sort by chromosome
  rows = sorted(csv.DictReader(sys.stdin, delimiter='\t'), key=lambda v: (int(v['Chr'][3:]) if (v['Chr'][3:]).isnumeric() else ord(v['Chr'][3:]), int(v['Start'])))
  logging.info('%i rows...', len(rows))
  
  for row in rows: 
    if chromosome is not None and row['Chr'] != chromosome:
      continue

    #Chr     Start   End     Gene    Mean    Median  Min     Max     Pct
    #chr17   0       77      none    311.5   297.0   138     483     0.038598
    data['max'].append(int(row['Max']))
    data['min'].append(int(row['Min']))
    data['mean'].append(float(row[colname]))

    # we started a new chromosome
    if row['Chr'] != last_chrom:
      chrom_offset = last_pos + 1
      if last_chrom is not None:
        ax.axvline(chrom_offset, alpha=0.5, color='#909090', lw=0.3)
        ax.text(chrom_offset, 0.99, last_chrom, color='#909090', ha='right', va='top', fontsize='x-small', rotation=90, transform=ax.get_xaxis_transform())
      last_chrom = row['Chr']
      logging.info('chromosome is now %s at %i', last_chrom, chrom_offset)
      #chromosomes[len(data['x'])] = row['Chr'] #last_pos # display on x-axis
      #chromosomes[pos] = row['Chr'] #last_pos # display on x-axis
    #elif highlight is not None and pos in highlight:
    #  ax.axvline(len(data['x']), alpha=0.5, color='#0000c0', lw=1)
    #  chromosomes[len(data['x'])] = pos
    #  #plt.text(len(data['x']) + 0.1, 0, row['Pos'], rotation=90)
    pos = chrom_offset + int(row['Start']) #+ (int(row['End']) - int(row['Start'])) / 2
    #logging.info('adding %s %i from start %s and offset %i', row['Chr'], pos, row['Start'], chrom_offset)
    data['x'].append(pos)
    last_pos = pos

  # normalise to total
  total_of_means = sum(data['mean'])
  data['mean'] = [x / total_of_means for x in data['mean']]
  #logging.info(data['mean'])

  # final chromosome
  #ax.axvline(chrom_offset, alpha=0.5, color='#909090', lw=0.3)
  ax.text(last_pos, 0.99, row['Chr'], color='#909090', ha='right', va='top', fontsize='x-small', rotation=90, transform=ax.get_xaxis_transform())

  if background is not None:
    logging.info('background %i files', len(background))
    background_means = []
    for n in background:
      means = []
      rows = sorted(csv.DictReader(open(n, 'rt'), delimiter='\t'), key=lambda v: (int(v['Chr'][3:]) if (v['Chr'][3:]).isnumeric() else ord(v['Chr'][3:]), int(v['Start'])))
      for row in rows:
        if chromosome is not None and row['Chr'] != chromosome:
          continue
        means.append(float(row[colname]))

      # normalise to total
      total_of_means = sum(means)
      means = [x / total_of_means for x in means]
      #logging.info(data['mean'])

      background_means.append(means)

    # plot mean, ci
    means = []
    sds = []
    lows = []
    highs = []
    for x in range(len(background_means[0])):
      vals = [background_means[n][x] for n in range(len(background_means))]
      sd = np.std(vals, ddof=1)
      mean = np.mean(vals)
      means.append(mean)
      sds.append(sd)
      lows.append(max(0, mean - 2 * sd))
      highs.append(mean + 2 * sd)

    if not normalise:
      ax.fill_between(data['x'], lows, highs, alpha=0.2, label='Background Range')
      ax.plot(data['x'], means, 'x-', label='Background Mean', linewidth=0.5, alpha=0.4, markersize=1, color='grey')
    #ax.plot(data['x'], data['min'], '-', label='Min', color='brown', linewidth=0.2, alpha=0.2)

  logging.info('writing to {}...'.format(target))
  #logging.info('chromosomes: %s', chromosomes)
  ax.grid(which='major', axis='y', linewidth=1)
  plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=True) # labels along the bottom edge are off

  if not normalise:
    ax.set_yscale('log', basey=2)
  if chromosome is None:
    ax.set_title(title)
  else:
    ax.set_title('{} - {}'.format(title, chromosome))
  plt.xticks(rotation=45, ha='right')
  #ax.plot(data['x'], data['max'], '-', label='Max', color='darkgreen', linewidth=0.2, alpha=0.2)
  #ax.fill_between(data['x'], data['min'], data['max'], alpha=0.2, label='Range')
  if normalise:
    #normalised = [(x[0] - x[1]) / x[2] for x in zip(data['mean'], means, sds)]
    normalised = [max(MAX_LOSS, math.log(max(1e-10, x[0]) / max(1e-10, x[1]), 2)) for x in zip(data['mean'], means, sds)]
    ax.plot(data['x'], normalised, 'o-', label='Log2ratio', linewidth=0.5, alpha=0.4, markersize=1, color='red')
    ax.set_ylabel(yaxis or 'Log2 of depth ratio to mean')

    # mark areas classified as gain or loss
    in_gain = in_loss = None
    for i in range(len(normalised)):
      if i < 2:
        continue
      if in_gain is not None and ((normalised[i] < GAIN and normalised[i-1] < GAIN) or i == len(normalised)-1): # two consecutive fails
        plt.axvspan(data['x'][in_gain], data['x'][i if i == len(normalised)-1 else i-2], facecolor=GAIN_COLOR, alpha=0.2)
        logging.info('exit gain from %i to %i', in_gain, i-2)
        in_gain = None
      elif in_gain is None and normalised[i] >= GAIN and (normalised[i-1] >= GAIN or normalised[i-2] >= GAIN): # two out of three positives
        in_gain = i-2 if normalised[i-2] >= GAIN else i-1
        logging.info('entered gain at %i', in_gain)        
      if in_loss is not None and ((normalised[i] > LOSS and normalised[i-1] > LOSS) or i == len(normalised)-1): # two consecutive fails
        plt.axvspan(data['x'][in_loss], data['x'][i if i == len(normalised)-1 else i-2], facecolor=LOSS_COLOR, alpha=0.2)
        logging.info('exit loss from %i to %i', in_loss, i-2)
        in_loss = None
      elif in_loss is None and normalised[i] <= LOSS and (normalised[i-1] <= LOSS or normalised[i-2] <= LOSS): # two out of three positives
        in_loss = i-2 if normalised[i-2] <= LOSS else i-1
        logging.info('entered loss at %i', in_loss)

  else:
    ax.plot(data['x'], data['mean'], 'o', label='Mean', linewidth=0.5, alpha=0.4, markersize=1, color='red')
    ax.set_ylabel(yaxis or 'Sequencing depth (reads)')
  #ax.plot(data['x'], data['min'], '-', label='Min', color='brown', linewidth=0.2, alpha=0.2)
  #ax.fill_between(x, y_est - y_err, y_est + y_err, alpha=0.2)

  ax.set_xlabel(xaxis or 'Genomic position')
  ax.legend()
  plt.tight_layout()
  fig.savefig(target, figsize=(width, height), dpi=dpi)
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Plot coverage')
  parser.add_argument('--title', required=False, default='Coverage', help='title for plot')
  parser.add_argument('--target', required=False, default='plot.png', help='target file')
  parser.add_argument('--highlight', required=False, nargs='*', help='positions to highlight')
  parser.add_argument('--xaxis', required=False, help='xaxis label')
  parser.add_argument('--yaxis', required=False, help='yaxis label')
  parser.add_argument('--chromosome', required=False, help='limit to this chromosome')
  parser.add_argument('--background', required=False, nargs='+', help='compare to these tumours')
  parser.add_argument('--normalise', action='store_true', help='normalise to mean and sd and background')
  parser.add_argument('--use_pct', action='store_true', help='use pct number of reads')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.target, args.title, args.highlight, xaxis=args.xaxis, yaxis=args.yaxis, chromosome=args.chromosome, background=args.background, normalise=args.normalise, pct=args.use_pct)
