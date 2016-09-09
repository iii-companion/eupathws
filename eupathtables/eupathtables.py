# -*- coding: utf-8 -*-
#
#  Author: Sascha Steinbiss <ss34@sanger.ac.uk>
#  Copyright (c) 2016 Genome Research Ltd
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
import six
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
    warnings.warn("GenomeTools Python bindings not found, " +
                  "stream interface not available")
    _gt_available = False


class WebServiceIterator(object):

    def _process_gene_model(self, table, gene):
        transcripts = {}
        for r in table['rows']:
            target_transcripts = None
            for f in r['fields']:
                if f['name'] == 'transcript_ids':
                    target_transcripts = f['value'].split(', ')
            assert(target_transcripts is not None)
            subfeat = {}
            for f in r['fields']:
                if f['name'] == 'type':
                    subfeat['Type'] = f['value'].lower()
                if f['name'] == 'sequence_id':
                    subfeat['Seqid'] = f['value']
                if f['name'] == 'gm_start':
                    subfeat['Start'] = int(f['value'])
                if f['name'] == 'gm_end':
                    subfeat['End'] = int(f['value'])
                if f['name'] == 'is_reversed':
                    subfeat['Strand'] = "+-"[int(f['value'])]
            for tt in target_transcripts:
                if not tt in transcripts:
                    transcripts[tt] = []
                transcripts[tt].append(subfeat)
        return transcripts

    def _process_pubmed(self, table, gene):
        pubmed_ids = []
        for r in table['rows']:
            for f in r['fields']:
                if f['name'] == 'pubmed_id':
                    pubmed_ids.append(f['value'])
        return pubmed_ids

    def _process_go(self, table, gene):
        goterms = []
        for r in table['rows']:
            target_transcripts = None
            for f in r['fields']:
                if f['name'] == 'transcript_ids':
                    target_transcripts = f['value'].split(', ')
            assert(target_transcripts is not None)
            goterm = {}
            # [GO ID] [Ontology]      [GO Term Name]  [Source]        [Evidence Code] [Is Not]
            for f in r['fields']:
                if f['name'] == 'ontology':
                    goterm['Ontology'] = f['value']
                if f['name'] == 'go_id':
                    goterm['GO ID'] = f['value']
                if f['name'] == 'go_term_name':
                    goterm['GO Term Name'] = f['value']
                if f['name'] == 'source':
                    goterm['Source'] = f['value']
                if f['name'] == 'evidence_code':
                    goterm['Evidence Code'] = f['value']
                if f['name'] == 'is_not':
                    goterm['Is Not'] = f['value']
            for tt in target_transcripts:
                if not tt in goterms:
                    goterms[tt] = []
                goterms[tt].append(goterm)
        return goterms

    def _json_to_gene(self, item):
        gene = {}
        gene['ID'] = item['id']

        for t in item['tables']:
            if t['name'] == 'Gene Model':
                gene['Gene Model'] = self._process_gene_model(t, gene)
            if t['name'] == 'PubMed':
                gene['PubMed'] = self._process_pubmed(t, gene)
            if t['name'] == 'GOTerms':
                gene['GO Terms'] = self._process_go(t, gene)

        return gene

    def _add_single_info(self, genes, item):
        fieldmap = {'sequence_id'  : 'seqid',
                    'gene_type'    : 'type',
                    'gene_product' : 'product',
                    'is_pseudo'    : 'pseudo'}
        if item['id'] not in genes:
            raise "gene not found with ID %s" % (item['id'])
        v = genes[item['id']]
        for f in item['fields']:
            if f['name'] in fieldmap:
                v[fieldmap[f['name']]] = f['value']
            elif f['name'] == 'gene_location_text':
                m = re.search('^([^:]+):([0-9,]+)\.\.([0-9,]+)\((.)\)', f['value'])
                if m:
                    v['start'] = int(m.group(2).replace(',',''))
                    v['stop'] = int(m.group(3).replace(',',''))
                    v['strand'] = m.group(4)
            else:
                v[f['name']] = f['value']

    def _get_cached_json(self, cache, url):
        j = None
        # get hostname part of base URL
        db = urllib.parse.urlparse(url).netloc
        if "ByTaxonGene.json" in url:
            cachefilename = db + "_" + self.organism.replace(" ", "_") + ".json"
        elif "ByTaxon.json" in url:
            cachefilename = db + "_" + self.organism.replace(" ", "_") + "_2.json"
        else:
            import md5
            cachefilename = db + "_" + self.organism.replace(" ", "_") + ("_%s.txt" % (md5.new(url).digest()))
        if not cache or (cache and not os.path.isfile(cachefilename)):
            res = requests.get(url, verify=True)
            if(res.ok):
                if cache:
                    with open(cachefilename, 'wb+') as f:
                        f.write(res.content)
                j = res.json()
            else:
                res.raise_for_status()
        else:
            with open(cachefilename, 'r') as f:
                j = json.loads(f.read())
        return j

    def __init__(self, baseurl, organism, cache=False):
        self.baseurl = baseurl
        self.organism = organism
        self.fields = ['annotated_go_function',
                       'gene_uniprot_id',
                       'gene_location_text',
                       'location_text',
                       'gene_type',
                       'gene_product',
                       'primary_key',
                       'ec_numbers',
                       'chromosome',
                       'gene_type',
                       'is_pseudo',
                       'source_id',
                       'sequence_id']
        self.tables = ['GeneModelDump',
                       'Notes',
                       'PubMed',
                       'GOTerms']
        url = ('{0}/webservices/GeneQuestions/GenesByTaxonGene.json?' +
               'organism={1}&o-tables={2}').format(self.baseurl, self.organism,
                                                   ','.join(self.tables))
        genes = {}

        # get gene models and transcripts
        taxongene_json = self._get_cached_json(cache, url)
        for v in taxongene_json['response']['recordset']['records']:
            gene = self._json_to_gene(v)
            if not gene['ID'] in genes:
                genes[gene['ID']] = gene
            else:
                raise 'duplicate gene ID: %s' % (gene['ID'])

        # get gene centric information
        url = ('{0}/webservices/GeneQuestions/GenesByTaxon.json?' +
               'organism={1}&o-fields={2}').format(self.baseurl, self.organism,
                                                   ','.join(self.fields))
        taxon_json = self._get_cached_json(cache, url)
        for v in taxon_json['response']['recordset']['records']:
            self._add_single_info(genes, v)
        self.genes = collections.deque(list(genes.values()))

    def next(self):
        return self.__next__()

    def __next__(self):
        if len(self.genes) > 0:
            return self.genes.popleft()
        else:
            raise StopIteration

if _gt_available:
    class TableInStream(gt.extended.CustomStream):

        def __init__(self, in_iterator, taxon_id=None):
            gt.extended.CustomStream.__init__(self)
            self.iterator = in_iterator

            self.outqueue = collections.deque()
            self.go_coll = GOCollection(taxon_id)
            self.uniprots = {}
            self.typemaps = {'gene': 'pseudogene',
                             'mRNA': 'pseudogenic_transcript',
                             'rRNA': 'pseudogenic_transcript',
                             'tRNA': 'pseudogenic_transcript',
                             'exon': 'pseudogenic_exon'}

        def get_transcript_length(self, model):
            transcript_length = 0
            for f in model:
                if int(f['Start']) <= int(f['End']) and f['Type'] == 'exon':
                    transcript_length += (int(f['End']) - int(f['Start']) + 1)
            return transcript_length

        def finaltype(self, v, mtype):
            if 'pseudo' in v and v['pseudo'] == 'Yes':
                if mtype in self.typemaps:
                    return self.typemaps[mtype]
            return mtype

        def make_noncoding(self, v, gtype, gene):
            if 'Gene Model' in v:
                for transcript_id, transcript_model in six.iteritems(v['Gene Model']):
                    transcript_length = self.get_transcript_length(transcript_model)
                    transcript = gt.extended.FeatureNode.create_new(v['seqid'], self.finaltype(v, gtype),
                                                                int(v['start']),
                                                                int(v['stop']),
                                                                v['strand'])
                    transcript.add_attribute("ID", transcript_id)
                    gene.add_child(transcript)
                    for f in transcript_model:
                        newfeat = gt.extended.FeatureNode.create_new(v['seqid'],
                                                                     self.finaltype(v, f['Type']),
                                                                     int(f['Start']),
                                                                     int(f['End']),
                                                                     v['strand'])
                        transcript.add_child(newfeat)


        def make_coding(self, v, gene):
            if 'Gene Model' in v:
                for transcript_id, transcript_model in six.iteritems(v['Gene Model']):
                    transcript_length = self.get_transcript_length(transcript_model)
                    transcript = gt.extended.FeatureNode.create_new(v['seqid'], self.finaltype(v, "mRNA"),
                                                                int(v['start']),
                                                                int(v['stop']),
                                                                v['strand'])
                    transcript.add_attribute("ID", transcript_id)
                    gene.add_child(transcript)

                    for f in transcript_model:
                        if int(f['Start']) <= int(f['End']):
                            newfeat = gt.extended.FeatureNode.create_new(f['Seqid'],
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
                                                  urllib.parse.quote(v['product'], safe=' '))
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
                    except StopIteration:
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
                    if 'gene_uniprot_id' in v:
                        for u in v['gene_uniprot_id'].split(','):
                            if len(u) > 0:
                                self.uniprots[u] = v['ID']

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
