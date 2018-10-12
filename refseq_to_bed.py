#!/usr/bin/env python
'''
  #bin    name    chrom   strand  txStart txEnd   cdsStart        cdsEnd  exonCount       exonStarts      exonEnds        score   name2   cdsStartStat    cdsEndStat      exonFrames
  0       NM_001308203.1  chr1    +       66999251        67216822        67000041        67208778        22      66999251,66999928,67091529,67098752,67105459,67108492,67109226,67136677,67137626,67138963,67142686,67145360,67154830,67155872,67160121,67184976,67194946,67199430,67205017,67206340,67206954,67208755,  66999355,67000051,67091593,67098777,67105516,67108547,67109402,67136702,67137678,67139049,67142779,67145435,67154958,67155999,67160187,67185088,67195102,67199563,67205220,67206405,67207119,67216822,  0       SGIP1   cmpl    cmpl    -1,0,1,2,0,0,1,0,1,2,1,1,1,0,1,1,2,2,0,2,1,1,
'''

import sys

EXCLUDE_UTR=True

for line in sys.stdin:
  if line.startswith('#'):
    continue

  fields = line.strip('\n').split('\t')
  chr = fields[2]
  gene = fields[12]
  cds_start = int(fields[6])
  cds_end = int(fields[7])
  for start, end in zip(fields[9].split(','), fields[10].split(',')):
    if start == '' or end == '':
      continue
    if EXCLUDE_UTR:
      start = max(int(start), cds_start)
      end = min(int(end), cds_end)
      if start < end:
        sys.stdout.write('{}\t{}\t{}\t{}\n'.format(chr, start, end, gene))
    else:
      sys.stdout.write('{}\t{}\t{}\t{}\n'.format(chr, start, end, gene))
