# EuPathTables
This package provides a Python interface for reading and converting EuPathDB 'gene information table' files as they are provided on the EuPathDB download site. The format in question is a structured text file in a custom format, containing most of the data available in the database in question. Here's an [example file](http://amoebadb.org/common/downloads/Current_Release/EhistolyticaHM1IMSS/txt/AmoebaDB-25_EhistolyticaHM1IMSSGene.txt).

## Usage
There are two ways of accessing the information in the file: via a Python generator returning one dict per gene, or via a GenomeTools input stream (which requires the [GenomeTools Python bindings](https://github.com/genometools/genometools/tree/master/gtpython)). This stream directly returns GenomeTools feature nodes for processing directly from the table without having to create GFF first.

Generator access:
```Python
#!/usr/bin/env python

import eupathtables

for g in eupathtables.iterate("../Ehist_PacBio/AmoebaDB-25_EhistolyticaHM1IMSSGene.txt"):
    print("%s\t%s:%s-%s" % (g['ID'], g['seqid'], g['start'], g['stop']))
```

Stream access:
```Python
#!/usr/bin/env python

import eupathtables
import gt

infile = "../Ehist_PacBio/AmoebaDB-25_EhistolyticaHM1IMSSGene.txt"
# we also create a GAF file with GO terms and products
gaf_out_file = "out.gaf"
# this is the taxon ID to use in the GAF file
taxon_id = 294381

table_in_stream = eupathtables.TableInStream(infile, gaf_out_file, taxon_id)
gff_out_stream = gt.extended.GFF3OutStream(table_in_stream)

fn = gff_out_stream.next_tree()
while fn:
    fn = gff_out_stream.next_tree()
```

We also provide a script for quick conversion to Companion-compatible GAF1 and GFF3:

```bash
  eupathtable_to_gff3 -g gaf.out -t 294381 ../Ehist_PacBio/AmoebaDB-25_EhistolyticaHM1IMSSGene.txt  > out.gff3
```

# Installation
Download/clone this repo from github, then:
```bash
python setup.py install
```

# Contact
ss34@sanger.ac.uk 
