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
import re
import requests
import os
import json
try:
    import exceptions
except:
    pass
import urllib
import sys
from eupathtables.go_collection import GOCollection
try:
    import gt
    _gt_available = True
except exceptions.ImportError:
    import warnings
    warnings.warn("GenomeTools Python bindings not found, stream interface not available")
    _gt_available = False

class WebServiceIterator(object):

    def _json_to_gene(self, item):
        gene = {}
        transcripts = {}
        gene['ID'] = item['id']
        for t in item['tables']:
            gene[t['name']] = {}
            for r in t['rows']:
                target_transcripts = None
                for f in r['fields']:
                    if f['name'] == 'transcript_ids':
                        target_transcripts = f['value'].split(', ')
                        print(target_transcripts)
                assert(target_transcripts is not None)
                for f in r['fields']:
                    subfeat = {}
                    if f['name'] == 'sequence_id':
                        subfeat['seqid'] = f['value']
                    if f['name'] == 'gm_start':
                        subfeat['start'] = int(f['value'])
                    if f['name'] == 'gm_end':
                        subfeat['end'] = int(f['value'])
                    if f['name'] == 'is_reversed':
                        subfeat['strand'] = "+-"[int(f['value'])]
                for tt in target_transcripts:
                    if not tt in transcripts:
                        transcripts[tt] = []
                    transcripts[tt].append(subfeat)
        gene.transcripts = transcripts

        return gene

    def __init__(self, database, organism, cache=False):
        self.db = database
        self.organism = organism
        self.cachefilename = self.db + "_" + organism.replace(" ", "_") + ".txt"
        self.fields = ['primary_key',
                       'chromosome',
                       'gene_type',
                       'is_pseudo',
                       'strand',
                       'source_id',
                       'sequence_id',
                       'five_prime_utr_length',
                       'three_prime_utr_length']
        self.tables = ['GeneModelDump',
                       'Notes',
                       'PubMed',
                       'GOTerms']
        url = ('http://{0}/webservices/GeneQuestions/GenesByTaxonGene.json?' +
               'organism={1}&o-tables={2}').format(self.db, self.organism,
                                                   ','.join(self.tables))
        if not cache or (cache and not os.path.isfile(self.cachefilename)):
            res = requests.get(url, verify=True)
            if(res.ok):
                if cache:
                    with open(cachefilename, 'wb+') as f:
                        f.write(res.content)
                j = res.json()
            else:
                res.raise_for_status()
        else:
            with open(self.cachefilename, 'r') as f:
                j = json.loads(f.read())
        for v in j['response']['recordset']['records']:
            print(self._json_to_gene(v))


class FlatFileIterator(object):

    def __init__(self, instream):
        self.instream = instream
        self.eof = False
        self.gene = dict()

    def __iter__(self):
        return self

    def next(self):
        return self.__next__()

    def __next__(self):
        self.gene = dict()
        while True:
            l = self.instream.readline()
            if l == '' or l is None:
                if not self.eof:
                    self.eof = True
                    return self.gene
                else:
                    raise exceptions.StopIteration
            m = re.match(r"Gene ID: (.*)$", l)
            if m and m.group(1):
                self.gene['ID'] = m.group(1)
                continue
            m = re.match(r"Gene Type: (.*)$", l)
            if m and m.group(1):
                self.gene['type'] = m.group(1)
                continue
            m = re.match(r"UniProt ID: (.*)$", l)
            if m and m.group(1) and m.group(1) != 'null':
                self.gene['uniprot_id'] = m.group(1)
                continue
            m = re.match(r"Gene Name: (.*)$", l)
            if m and m.group(1) and m.group(1) != 'null':
                self.gene['name'] = m.group(1)
                continue
            m = re.match(r"Is Pseudo: (.*)$", l)
            if m and m.group(1):
                if m.group(1) == 'Yes':
                    self.gene['pseudo'] = True
                else:
                    self.gene['pseudo'] = False
                continue
            m = re.match(r"Annotated ([35]). UTR length: ([0-9]+)", l)
            if m and m.group(1) and m.group(2):
                self.gene["utr_%s" % m.group(1)] = int(m.group(2))
                continue
            m = re.match(r"Product Description: (.*)$", l)
            if m and m.group(1):
                self.gene['product'] = m.group(1)
                continue
            m = re.match(r"Genomic Location: (?P<seqid>[^ ]+): " +
                         "(?P<start>[0-9,]+) - (?P<stop>[0-9,]+) " +
                         "\((?P<strand>.)\)", l)
            if m:
                for k, v in m.groupdict().iteritems():
                    self.gene[k] = v.replace(',', '')
                continue
            m = re.match(r"TABLE: (.*)$", l)
            if m and m.group(1):
                thistable = []
                tablename = m.group(1)
                l = self.instream.readline()
                listidx = re.findall("\[([^\]]+)\]", l)
                l = self.instream.readline()
                while l != '\n' and l != '':
                    i = 0
                    thisline = dict()
                    for v in re.split('\t', l):
                        thisline[listidx[i]] = v.rstrip()
                        i = i + 1
                    thistable.append(thisline)
                    l = self.instream.readline()
                self.gene[tablename] = thistable
                continue
            if re.match('------------------', l):
                return self.gene


if _gt_available:
    class TableInStream(gt.extended.CustomStream):

        def __init__(self, in_iterator, taxon_id=None):
            gt.extended.CustomStream.__init__(self)
            if not it:
                self.iterator = FlatFileIterator(instream)
            else:
                self.iterator = it
            self.outqueue = collections.deque()
            self.go_coll = GOCollection(taxon_id)
            self.uniprots = {}
            self.typemaps = {'gene': 'pseudogene',
                             'mRNA': 'pseudogenic_transcript',
                             'rRNA': 'pseudogenic_transcript',
                             'tRNA': 'pseudogenic_transcript',
                             'exon': 'pseudogenic_exon'}

        def get_transcript_length(self, v):
            transcript_length = 0
            for f in v['Gene Model']:
                if int(f['Start']) <= int(f['End']) and f['Type'] == 'exon':
                    transcript_length += (int(f['End']) - int(f['Start']) + 1)
            return transcript_length

        def finaltype(self, v, mtype):
            if 'pseudo' in v and v['pseudo']:
                if mtype in self.typemaps:
                    return self.typemaps[mtype]
            return mtype

        def make_noncoding(self, v, gtype, gene):
            rna = gt.extended.FeatureNode.create_new(v['seqid'], self.finaltype(v, gtype),
                                                     int(v['start']), int(v['stop']),
                                                     v['strand'])
            rna.add_attribute("ID", v['ID'] + ":" + gtype)
            gene.add_child(rna)
            if 'Gene Model' in v:
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

        def make_coding(self, v, gene):
            transcript = gt.extended.FeatureNode.create_new(v['seqid'], self.finaltype(v, "mRNA"),
                                                            int(v['start']),
                                                            int(v['stop']),
                                                            v['strand'])
            transcript.add_attribute("ID", v['ID'] + ".1")
            gene.add_child(transcript)
            if 'Gene Model' in v:
                newend = 1
                left_i = 0
                left_type = 'five_prime_UTR'
                left_key = 'utr_5'
                right_type = 'three_prime_UTR'
                right_key = 'utr_3'
                left_c = None
                right_c = None

                if 'utr_5' in v or 'utr_3' in v:
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

                transcript_length = self.get_transcript_length(v)

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
            if 'product' in v and v['product'] != "null":
                polypeptide.add_attribute("product",
                                          'term%3D' +
                                          urllib.quote(v['product'], safe=' '))
            self.outqueue.append(polypeptide)

            # register GO terms
            if self.go_coll:
                self.go_coll.add_item(v)

        def next(self):
            while True:
                #try:
                    if len(self.outqueue) > 0:
                        return self.outqueue.popleft()

                    try:
                        v = self.iterator.next()
                    except exceptions.StopIteration:
                        return None
                    if not v:
                        return None

                    gene = gt.extended.FeatureNode.create_new(v['seqid'], self.finaltype(v, "gene"),
                                                              int(v['start']),
                                                              int(v['stop']), v['strand'])
                    gene.add_attribute("ID", v['ID'])
                    if 'name' in v:
                        gene.add_attribute("Name", v['name'])

                    # track UniProtID -> gene mappings
                    # XXX handle multiple transcripts per product once supported by
                    # EuPathDB!
                    if 'uniprot_id' in v:
                        self.uniprots[v['uniprot_id']] = v['ID']

                    # non-coding RNA
                    m = re.match(r"(.RNA) encoding", v['type'])
                    if m:
                        self.make_noncoding(v, m.group(1), gene)

                    # protein coding gene
                    if v['type'] == 'protein coding':
                        self.make_coding(v, gene)

                    break

                #except:
                #    if v:
                #        sys.stderr.write("error creating feature for %s\n" % v['ID'])
                #    else:
                #        sys.stderr.write("error creating feature , no ID yet\n")
                #    continue

            return gene
