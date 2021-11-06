# THIS IS NOT WORKING YET - DO NOT USE

# Pylign

Pylign is an experimental package to simplify realignment to a reference for viral sequences, using only Python packages. It owes almost everything to [sam_2_fasta](https://github.com/cov-ert/datafunk/blob/master/datafunk/sam_2_fasta.py) by Ben Jackson at the Univesity of Edinburgh. (Sam2fasta has since been ported to [gofasta](https://github.com/cov-ert/gofasta)). Sam_2_fasta is typically run after `minimap2` alignment. Pylign incorporates the alignment process, using [mappy](https://pypi.org/project/mappy/).

