
# Papillon
![Papillon](assets/papillon.png)

## Installation
This software requires Python3.

```
python -m venv papillon-env
source ./papillon-env/bin/activate
pip install -r requirements.txt
```

## Usage

### Generate a bed file containing regions of interest

For example:
```
refseq_to_bed.py --coding --genes MSH2 MSH6 MLH1 PMS2 --transcripts NM_000251 NM_000179 NM_000249 NM_000535 --refseq  /data/projects/punim0567/data/public_datasets/refseq.ucsc.hg19.180829.gz
```

* the coding flag excludes UTR regions
* the genes argument (optional) chooses the genes to include
* the transcripts argument (optional) chooses the transcripts to include
* the output is a bed file

### Find coverage across regions of interest and generate heatmaps
Generate plots and a table of coverage across the exons in the specified genes:

For example:
```
papillon.py \
    --max_coverage 100 \
    --stat median \
    --bed refseq.bed \
    --bams *.bam \
    --genes MSH2 MSH6 MLH1 PMS2 \
    --plot gene_coverage.tumour.median \
    --gene_level gene_coverage.tumour.median.tsv \
    --padding 20 \
    --exon_plots \
    --exon_number \
    --exon_number_reverse PMS2 \
    --sample_name_end '-' \
    > out/papillon.tumour.median.tsv &
```

* max_coverage limits the maximum value shown on the plot
* stat is the statistic to report on the plot
* bed is a bed file containing regions to examine. If --genes is specified it will filter on the fourth column.
* plot is the main plot name prefix
* gene_level (optional) generates a tsv with gene level statistics
* padding extends the regions on either side
* exon_plots generates plots for each gene with individual regions
* exon_number uses incremental region number on the x-axis
* exon_number_reverse for genes on the reverse strand
* sample_name_end uses as the sample name everything in the bam up to this character

### Summarise into gene

```
python summary.py --summary overall.tsv < heat.tsv > summary.tsv
```

### Find potential deletions or duplications
```
python outliers.py < heat.tsv > outliers.tsv
```

## Authors
* Peter Georgeson
* Romy Walker
