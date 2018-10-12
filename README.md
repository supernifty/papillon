
# Papillon

## Installation
This software requires Python3.

```
python -m venv papillon-env
source ./papillon-env/bin/activate
pip install -r requirements.txt
```

## Usage

```
python papillon.py --bed refseq.bed --bams *.bam --genes MSH2 MSH3 MSH6 MLH1 --plot batch1 > heat1.tsv
```

## Authors
* Peter Georgeson
* Romy Walker
