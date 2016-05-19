# EuPathTables [![Build Status](https://api.travis-ci.org/sanger-pathogens/eupathtables.svg)](https://travis-ci.org/sanger-pathogens/eupathtables)
This package provides a Python interface for reading and converting EuPathDB 'gene information table' files as they are provided on the EuPathDB download site. The format in question is a structured text file in a custom format, containing most of the data available in the database in question. Here's an [example file](http://fungidb.org/common/downloads/release-28/Aniger_ATCC1015/txt/FungiDB-28_Aniger_ATCC1015Gene.txt).

EuPathTables also recognises UTRs and pseudogenes and provides this information in appropriate fields/types.

## Usage
There are two ways of accessing the information in the file: via a Python iterator returning one dict per gene, or via a GenomeTools input stream (which requires the [GenomeTools Python bindings](https://github.com/genometools/genometools/tree/master/gtpython)). This stream directly returns GenomeTools feature nodes for processing directly from the table without having to create GFF first.

Generator access:
```Python
#!/usr/bin/env python

import eupathtables

for g in eupathtables.FlatFileIterator(open('FungiDB-28_Aniger_ATCC1015Gene.txt')):
    print("%s\t%s:%s-%s" % (g['ID'], g['seqid'], g['start'], g['stop']))
```

Stream access:
```Python
#!/usr/bin/env python

import eupathtables
import gt

infile = "FungiDB-28_Aniger_ATCC1015Gene.txt"
# we also create a GAF file with GO terms and products
gaf_out_file = "out.gaf"
# this is the taxon ID to use in the GAF file
taxon_id = 294381

table_in_stream = eupathtables.TableInStream(open(infile), taxon_id)
gff_out_stream = gt.extended.GFF3OutStream(table_in_stream)

fn = gff_out_stream.next_tree()
while fn:
    fn = gff_out_stream.next_tree()

# write GO terms out to GAF 1.0 file
table_in_stream.go_coll.to_gafv1(open(gaf_out_file, "w+"))
```

We also provide a script for quick conversion to Companion-compatible GAF1 and GFF3:

```bash
  eupathtable_to_gff3 -g gaf.out -t 294381 FungiDB-28_Aniger_ATCC1015Gene.txt  > out.gff3
```

# Installation
Download/clone this repo from github, then:
```bash
python setup.py install
```

# Contact
ss34@sanger.ac.uk
