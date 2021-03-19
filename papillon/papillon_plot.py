#!/usr/bin/env python
'''
  measure coverage across regions and bams
'''

import argparse
import csv
import collections
import logging
import sys

import intervaltree
import numpy as np
import pysam
import matplotlib as mpl
from matplotlib.ticker import FuncFormatter
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

def fmt_container(m):
  def x_fmt(x, pos):
    if x in m:
      return m[x]
    else:
      return ''
  return x_fmt

def main(target, title, highlight, dpi=300, width=12, height=8):
  logging.info('starting...')

  last_pos = None
  data = {'x': [], 'max': [], 'sdu': [], 'mean': [], 'sdl': [], 'min': []} # max, sdu, mean, sdl, min
  exons = {}
  fig, ax = plt.subplots()
  ax.xaxis.set_major_formatter(FuncFormatter(fmt_container(exons)))
  ax.xaxis.set_minor_formatter(FuncFormatter(fmt_container(exons)))

  for row in csv.DictReader(sys.stdin, delimiter='\t'): 
    # Chr     Pos     n       Mean    Median  Min     Max     SD
    #chr1    45794923        3093    196.1   184     9       1667    158.6
    data['max'].append(int(row['Max']))
    data['min'].append(int(row['Min']))
    data['mean'].append(float(row['Mean']))
    data['sdu'].append(float(row['pct97_5']))
    data['sdl'].append(float(row['pct2_5']))
    if last_pos is not None and int(row['Pos']) - last_pos > 1:
      ax.axvline(len(data['x']), alpha=0.1, color='#909090', lw=1)
      #exons[len(data['x'])] = last_pos # display on x-axis
    elif row['Pos'] in highlight:
      ax.axvline(len(data['x']), alpha=0.5, color='#0000c0', lw=1)
      exons[len(data['x'])] = int(row['Pos'])
      #plt.text(len(data['x']) + 0.1, 0, row['Pos'], rotation=90)
    data['x'].append(row['Pos'])
    last_pos = int(row['Pos'])


  logging.info('writing to {}...'.format(target))
  ax.grid(which='major', axis='y', linewidth=1)
  plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=True) # labels along the bottom edge are off
  ax.set_yscale('log')
  ax.set_title(title)
  plt.xticks(rotation=45, ha='right')
  ax.plot(data['x'], data['max'], '-', label='Max', color='darkgreen')
  ax.fill_between(data['x'], data['sdl'], data['sdu'], alpha=0.2, label='95%')
  ax.plot(data['x'], data['mean'], '-', label='Mean')
  ax.plot(data['x'], data['min'], '-', label='Min', color='brown')
  #ax.fill_between(x, y_est - y_err, y_est + y_err, alpha=0.2)

  ax.set_xlabel('Position')
  ax.set_ylabel('Depth')
  ax.legend()
  plt.tight_layout()
  fig.savefig(target, figsize=(width, height), dpi=dpi)
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Plot coverage')
  parser.add_argument('--title', required=False, default='Coverage', help='title for plot')
  parser.add_argument('--target', required=False, default='plot.png', help='target file')
  parser.add_argument('--highlight', required=False, nargs='*', help='positions to highlight')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.target, args.title, args.highlight)
