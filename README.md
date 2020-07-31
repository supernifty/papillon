
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

### Find coverage across regions of interest
Generate plots and a table of coverage across the exons in the specified genes:
```
python papillon.py --bed refseq.bed --bams *.bam --genes MSH2 MSH3 MSH6 MLH1 --plot batch1 > heat.tsv
```

* refseq.bed is a bed file containing exons to examine. If --genes is specified it will filter on the fourth column.

Summarise into gene
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
