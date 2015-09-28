#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Author: Sascha Steinbiss <ss34@sanger.ac.uk>
#  Copyright (c) 2015 Genome Research Ltd
#
#  Permission to use, copy, modify, and distribute this software for any
#  purpose with or without fee is hereby granted, provided that the above
#  copyright notice and this permission notice appear in all copies.
#
#  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
#  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
#  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
#  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
#  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
#  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
#  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

import collections
import datetime
import gt
import pprint
import re
import urllib
import sys

today = datetime.date.today()

def iterate(table_file):
    myfile = open(table_file)
    gene = dict()
    while True:
        l = myfile.readline()
        if l == '':
            yield None
        m = re.match(r"Gene ID: (.*)$", l)
        if m and m.group(1):
            gene['ID'] = m.group(1)
            continue
        m = re.match(r"Gene Type: (.*)$", l)
        if m and m.group(1):
            gene['type'] = m.group(1)
            continue
        m = re.match(r"Gene Name: (.*)$", l)
        if m and m.group(1):
            gene['name'] = m.group(1)
            continue
        m = re.match(r"Product Description: (.*)$", l)
        if m and m.group(1):
            gene['product'] = m.group(1)
            continue
        m = re.match(r"Genomic Location: (?P<seqid>[^ ]+): " +
                      "(?P<start>[0-9,]+) - (?P<stop>[0-9,]+) " +
                      "\((?P<strand>.)\)", l)
        if m:
            for k, v in m.groupdict().iteritems():
                gene[k] = v.replace(',', '')
            continue
        m = re.match(r"TABLE: (.*)$", l)
        if m and m.group(1):
            thistable = []
            tablename = m.group(1)
            l = myfile.readline()
            listidx = re.findall("\[([^\]]+)\]", l)
            l = myfile.readline()
            while l != '\n':
                i = 0
                thisline = dict()
                for v in re.split('\t', l):
                    thisline[listidx[i]] = v.rstrip()
                    i = i + 1
                thistable.append(thisline)
                l = myfile.readline()
            gene[tablename] = thistable
            continue
        if re.match('------------------', l):
            yield gene
            gene = dict()

class TableInStream(gt.extended.CustomStream):

    def __init__(self, table_file, gaf_file=None, taxon_id=None):
        gt.extended.CustomStream.__init__(self)
        self.iterator = iterate(table_file)
        self.outqueue = collections.deque()
        self.aspects = {'Biological Process': 'P',
                        'Molecular Function': 'F',
                        'Cellular Component': 'C'}
        if gaf_file:
          self.gaf_file = open(gaf_file, 'w+')
          self.gaf_file.write("!gaf-version: 1.0\n")
          self.taxon_id = taxon_id

    def aspect2oneletter(self, aspect):
        return self.aspects[aspect]

    def next(self):
        if len(self.outqueue) > 0:
            return self.outqueue.popleft()
        v = self.iterator.next()

        if not v:
            return None
        # make gene
        gene = gt.extended.FeatureNode.create_new(v['seqid'], "gene",
                                                  int(v['start']),
                                                  int(v['stop']), v['strand'])
        gene.add_attribute("ID", v['ID'])
        if v['name'] != 'null':
            gene.add_attribute("Name", v['name'])
        if v['type'] == 'protein coding':
            # make transcript
            transcript = gt.extended.FeatureNode.create_new(v['seqid'], "mRNA",
                                                            int(v['start']),
                                                            int(v['stop']),
                                                            v['strand'])
            transcript.add_attribute("ID", v['ID'] + ".1")
            gene.add_child(transcript)
            # make introns/exons/...
            if v['Gene Model']:
                for f in v['Gene Model']:
                    newfeat = gt.extended.FeatureNode.create_new(v['seqid'],
                                                                 f['Type'],
                                                                 int(f['Start']),
                                                                 int(f['End']),
                                                                 v['strand'])
                    transcript.add_child(newfeat)
                    if f['Type'] == 'exon':
                        newfeat = gt.extended.FeatureNode.create_new(v['seqid'],
                                                                'CDS',
                                                                int(f['Start']),
                                                                int(f['End']),
                                                                v['strand'])
                        transcript.add_child(newfeat)

            # make polypeptide
            polypeptide = gt.extended.FeatureNode.create_new(v['seqid'],
                                                             "polypeptide",
                                                             int(v['start']),
                                                             int(v['stop']),
                                                             v['strand'])
            polypeptide.add_attribute("ID",
                                      transcript.get_attribute("ID") + ":pep")
            polypeptide.add_attribute("Derives_from",
                                      transcript.get_attribute("ID"))
            if v['product'] != "null":
                polypeptide.add_attribute("product",
                                          'term%3D' +
                                          urllib.quote(v['product'], safe=' '))
            self.outqueue.append(polypeptide)

            if hasattr(self, 'gaf_file') and v['GO Terms']:
                terms_used = {}
                for go in v['GO Terms']:
                    aspect = None
                    try:
                        aspect = self.aspect2oneletter(go['Ontology'])
                    except:
                        sys.stderr.write("error trying to get aspect for: "
                                         + str(go) + "\n")
                        continue
                    if aspect and not go['GO ID'] in terms_used:
                        self.gaf_file.write("EuPathDB\t" + str(v['ID']) + "\t"
                                        + str(v['ID']) + "\t\t" + go['GO ID']
                                        + "\tGO_REF:0000001\t"
                                        + go['Evidence Code'] + "\t\t"
                                        + aspect + "\t" + v['product'] + "\t\t"
                                        + "gene\ttaxon:" + str(self.taxon_id)
                                        + "\t"
                                        + today.strftime('%Y%m%d')
                                        + "\tEuPathDB\n")
                        terms_used[go['GO ID']] = True

        return gene
