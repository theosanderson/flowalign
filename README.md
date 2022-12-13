# FlowAlign

FlowAlign is a package to simplify realignment to a reference for viral sequences, using only Python packages. (Under the hood it is using minimap2, so it's fast.)

## Installation

```bash
pip install flowalign
```
## Usage (command-line)
```bash
flowalign sequences.fasta --reference wuhCor1.fa --output aligned.fa
```

If you omit the output parameter, output will go to STDOUT.

## Usage (Python):
First we download the reference, and some unaligned sequences to align to it:
```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/bigZips/wuhCor1.fa.gz && gunzip wuhCor1.fa.gz
wget https://data.nextstrain.org/files/ncov/open/global/sequences.fasta.xz &&  xz --decompress sequences.fasta.xz
```

Then we write a simple Python script, creating an iterator called `aligned` that will yield reference-aligned versions of these unaligned sequences:
```py
import flowalign
aligned = flowalign.yield_aligned(input="sequences.fasta", reference= "wuhCor1.fa")
for name, aligned_sequence in aligned:
    print(">"+name)
    print(aligned_sequence)
```

yield_aligned can also take a stream, e.g.:
```py
aligned = flowalign.yield_aligned(input=open("sequences.fasta","rt), reference= "wuhCor1.fa")
```

Under the hood, mappy is in some sense calling minimap2 with `--secondary=no --sam-hit-only --score-N=0 -x asm5`.

Note that the multiprocessing implementation is fairly hacky which may cause issues. It is mostly expected that `yield_aligned` will be only called once at any one time.




## Acknowledgements

This package is developed by me, but it owes almost everything to [sam_2_fasta](https://github.com/cov-ert/datafunk/blob/master/datafunk/sam_2_fasta.py) by [Ben Jackson](https://github.com/benjamincjackson) at the University of Edinburgh -- [major functions](https://github.com/theosanderson/flowalign/blob/main/src/flowalign/functions_based_on_sam_2_fasta.py) are taken from that codebase. (Ben has also ported this code to [gofasta](https://github.com/cov-ert/gofasta)). Sam2fasta is typically run on a SAM file from [minimap2](https://github.com/lh3/minimap2). flowalign incorporates the alignment process using [mappy](https://pypi.org/project/mappy/), the Python-bindings for minimap2. The idea is that one doesn't need any dependencies except Python packages (mappy supplies its own minimap2) to get aligned sequences.

