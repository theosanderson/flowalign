# Pylign

Pylign is an experimental package to simplify realignment to a reference for viral sequences, using only Python packages. This package is developed by me, but it owes almost everything to [sam_2_fasta](https://github.com/cov-ert/datafunk/blob/master/datafunk/sam_2_fasta.py) by [Ben Jackson](https://github.com/benjamincjackson) at the University of Edinburgh. (Ben has also ported this code to [gofasta](https://github.com/cov-ert/gofasta)). Sam2fasta is typically run on a SAM file from [minimap2](https://github.com/lh3/minimap2). Pylign incorporates the alignment process using [mappy](https://pypi.org/project/mappy/), the Python-bindings for minimap2. The idea is that one doesn't need any dependencies except Python packages (mappy supplies its own minimap2) to get aligned sequences.

## Installation

```bash
pip install pylign
```


## Usage (Python):
First we download the reference, and some unaligned sequences to align to it:
```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/bigZips/wuhCor1.fa.gz && gunzip wuhCor1.fa.gz
wget https://data.nextstrain.org/files/ncov/open/global/sequences.fasta.xz &&  xz --decompress sequences.fasta.xz
```

Then we write a simple Python script, creating an iterator called `aligned` that will yield reference-aligned versions of these unaligned sequences:
```py
import pylign
aligned = pylign.yield_aligned(input="sequences.fasta", reference= "wuhCor1.fa")
for name, aligned_sequence in aligned:
    print(">"+name)
    print(aligned_sequence)
```

yield_aligned can also take a stream, e.g.:
```py
aligned = pylign.yield_aligned(input=open("sequences.fasta","rt), reference= "wuhCor1.fa")
```

Under the hood, mappy is in some sense calling minimap2 mode with `--secondary=no --sam-hit-only --score-N=0 -x asm20`.

Note that the multiprocessing implementation is fairly hacky which may cause issues. Expected use is that this function will be only called once, and certainly only once at any particular time.

## Usage (command-line)
```bash
pylign sequences.fasta --reference wuhCor1.fa --output aligned.fa
```

If you omit the output, output will go to STDOUT.

