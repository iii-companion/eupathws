#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Author: Sascha Steinbiss <ss34@sanger.ac.uk>
#  Copyright (c) 2015-2016 Genome Research Ltd
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
        if l == '' or l == None:
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
        m = re.match(r"Is Pseudo: (.*)$", l)
        if m and m.group(1):
            if m.group(1) == 'Yes':
                gene['pseudo'] = True
            else:
                gene['pseudo'] = False
            continue
        m = re.match(r"Annotated ([35]). UTR length: ([0-9]+)", l)
        if m and m.group(1) and m.group(2):
            gene["utr_%s" % m.group(1)] = int(m.group(2))
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
        self.typemaps = {'gene': 'pseudogene',
                         'mRNA': 'pseudogenic_transcript',
                         'rRNA': 'pseudogenic_transcript',
                         'tRNA': 'pseudogenic_transcript',
                         'exon': 'pseudogenic_exon'}
        if gaf_file:
          self.gaf_file = open(gaf_file, 'w+')
          self.gaf_file.write("!gaf-version: 1.0\n")
          self.taxon_id = taxon_id

    def aspect2oneletter(self, aspect):
        return self.aspects[aspect]

    def finaltype(self, v, mtype):
        if v['pseudo']:
            if mtype in self.typemaps:
                return self.typemaps[mtype]
        return mtype

    def next(self):
        # make gene
        v = None
        while True:
            #try:
                if len(self.outqueue) > 0:
                    return self.outqueue.popleft()
                v = self.iterator.next()

                if not v:
                    return None

                gene = gt.extended.FeatureNode.create_new(v['seqid'], self.finaltype(v, "gene"),
                                                          int(v['start']),
                                                          int(v['stop']), v['strand'])
                gene.add_attribute("ID", v['ID'])
                if v['name'] != 'null':
                    gene.add_attribute("Name", v['name'])

                # non-coding RNA
                m = re.match(r"(.RNA) encoding", v['type'])
                if m:
                    # make RNA
                    rna = gt.extended.FeatureNode.create_new(v['seqid'], self.finaltype(v, m.group(1)),
                                                             int(v['start']), int(v['stop']),
                                                             v['strand'])
                    rna.add_attribute("ID", v['ID'] + ":" + m.group(1))
                    gene.add_child(rna)
                    # make introns/exons/...
                    if v['Gene Model']:
                        for f in v['Gene Model']:
                            if int(f['Start']) <= int(f['End']):
                                newfeat = gt.extended.FeatureNode.create_new(v['seqid'],
                                                                             self.finaltype(v, f['Type']),
                                                                             int(f['Start']),
                                                                             int(f['End']),
                                                                             v['strand'])
                                rna.add_child(newfeat)
                            else:
                                sys.stderr.write("invalid feature range, skipping "
                                                 + str(f['Type']) + " " + v['ID']
                                                 + "\n")

                # protein coding gene
                if v['type'] == 'protein coding':
                    # make transcript
                    transcript = gt.extended.FeatureNode.create_new(v['seqid'], self.finaltype(v, "mRNA"),
                                                                    int(v['start']),
                                                                    int(v['stop']),
                                                                    v['strand'])
                    transcript.add_attribute("ID", v['ID'] + ".1")
                    gene.add_child(transcript)

                    # make introns/exons/CDS...
                    if v['Gene Model']:
                        newend = 1
                        left_i = 0
                        left_type = 'five_prime_UTR'
                        left_key = 'utr_5'
                        right_type = 'three_prime_UTR'
                        right_key = 'utr_3'
                        left_c = None
                        right_c = None

                        if 'utr_5' in v or 'utr_3' in v:
                            #if 'utr_3' in v:
                            #    print(">>>> %s (%s) has 3 prime UTR %s" % (v['ID'], v['strand'], str(transcript_length - v['utr_3'])))
                            if v['strand'] == '-':
                                if 'utr_5' in v:
                                    right_c = int(v['utr_5'])
                                if 'utr_3' in v:
                                    left_c = int(v['utr_3'])
                                left_type = 'three_prime_UTR'
                                right_type = 'five_prime_UTR'
                                left_key = 'utr_3'
                                right_key = 'utr_5'
                            else:
                                if 'utr_5' in v:
                                    left_c = int(v['utr_5'])
                                if 'utr_3' in v:
                                    right_c = int(v['utr_3'])

                        if v['strand'] == '-':
                            v['Gene Model'].reverse()

                        transcript_length = 0
                        for f in v['Gene Model']:
                            if int(f['Start']) <= int(f['End']) and f['Type'] == 'exon':
                                transcript_length += (int(f['End']) - int(f['Start']) + 1)

                        for f in v['Gene Model']:
                            if int(f['Start']) <= int(f['End']):
                                newfeat = gt.extended.FeatureNode.create_new(v['seqid'],
                                                                             self.finaltype(v, f['Type']),
                                                                             int(f['Start']),
                                                                             int(f['End']),
                                                                             v['strand'])
                                transcript.add_child(newfeat)

                                # decide whether to make CDS
                                if f['Type'] == 'exon':
                                    scoord = int(f['Start'])
                                    ecoord = int(f['End'])
                                    newend = left_i + (int(f['End']) - int(f['Start']) + 1)
                                    #print(">>>> %s: newend at %d" % (v['ID'], newend))

                                    if left_key in v:
                                        if newend > left_c:
                                            if left_i < left_c:
                                                newfeat = gt.extended.FeatureNode.create_new(v['seqid'],
                                                                 self.finaltype(v, left_type),
                                                                 scoord,
                                                                 int(f['Start']) + (left_c - left_i) - 1,
                                                                 v['strand'])
                                                transcript.add_child(newfeat)
                                                scoord = int(f['Start']) + (left_c - left_i)
                                                #print("left_i %d left_c %d" % (left_i, left_c))
                                            else:
                                                scoord = int(f['Start'])
                                        else:
                                            newfeat = gt.extended.FeatureNode.create_new(v['seqid'],
                                                                 self.finaltype(v, left_type),
                                                                 scoord, ecoord,
                                                                 v['strand'])
                                            transcript.add_child(newfeat)
                                            scoord = None

                                    if right_key in v:
                                        if newend > (transcript_length - right_c):
                                            #print(">>>> %s: crossed 3 prime UTR line at %d" % (v['ID'], newend))
                                            if left_i < (transcript_length - right_c):
                                                newfeat = gt.extended.FeatureNode.create_new(v['seqid'],
                                                                 self.finaltype(v, right_type),
                                                                 int(f['Start']) + ((transcript_length - right_c) - left_i),
                                                                 ecoord,
                                                                 v['strand'])
                                                transcript.add_child(newfeat)
                                                ecoord = int(f['Start']) + ((transcript_length - right_c) - left_i - 1)
                                            else:
                                                newfeat = gt.extended.FeatureNode.create_new(v['seqid'],
                                                                 self.finaltype(v, right_type),
                                                                 scoord, ecoord,
                                                                 v['strand'])
                                                transcript.add_child(newfeat)
                                                ecoord = None

                                    if ecoord and scoord:
                                        newfeat = gt.extended.FeatureNode.create_new(v['seqid'],
                                                                     self.finaltype(v, 'CDS'),
                                                                     scoord, ecoord,
                                                                     v['strand'])
                                        transcript.add_child(newfeat)
                                    left_i = newend
                                #print(">>> %s left_i %d" % (v['ID'],left_i))

                            else:
                                sys.stderr.write("invalid feature range, skipping "
                                                 + str(f['Type']) + " " + v['ID']
                                                 + "\n")

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

                break
            # except:
            #     if v:
            #         sys.stderr.write("error creating feature for %s\n" % v['ID'])
            #     else:
            #         sys.stderr.write("error creating feature , no ID yet\n")
            #     continue

        return gene