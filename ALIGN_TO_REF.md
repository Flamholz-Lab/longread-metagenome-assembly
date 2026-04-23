# Notes for aligning your assembled genome to a reference

## There are a few tools for this. Here are the ones I've tested:
ec_reference.fasta is my reference genome and ec_meta_assemble.fa is what I assembled from ONT sequencing

### 1. Mummer4
Mummer can be installed via conda. I used:

```
nucmer --prefix=ecoli ec_reference.fasta ec_meta_assemble.fa
dnadiff -d ecoli.delta
```

and got a nice out.report that shows alignment stats, number of SNPs, insertions, etc. This is a good summary.

Then to look at potential insertions in the assembled genome relative to the reference, I used:

```
show-diff -q ecoli.delta > differences.txt

python parse_dnadiff.py \
    -d differences.txt \
    -a ec_meta_assemble.fa \
    -o candidates.fna
```

And used Prodigal (conda-installable) to translate the results:

```
prodigal -i candidates.fna \
    -a proteins.faa \
    -o genes.gff \
    -f gff \
    -p meta
```

Then I blasted the longest results to look at ponetial inserted genes in the assembled genome relative to the NCBI reference.


### 2. QUAST
QUAST is specifically designed for assemblies vs. a reference. I had a little trouble installing, here's how I did it:

```
conda create -n quast
conda activate quest
conda install conda-forge::zlib
conda install pip 
pip install quast
pip install setuptools
```

And then running:

```
quast.py ./dastool_output/SemiBin_1.fa -r ec_reference.fasta -o quast_output
```

This provides a nice .html visual of your contigs vs. the reference genome, as well as GC content, but doesn't provide the same stats or a full alignment like mummer.
