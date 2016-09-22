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
from eupathtables.login import get_session, parse_login
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import six
import json
import urllib
import sys
from eupathtables.go_collection import GOCollection
try:
    import gt
    _gt_available = True
except ImportError:
    import warnings
    warnings.warn("GenomeTools Python bindings not found, " +
                  "stream interface not available")
    _gt_available = False
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class WebServiceIterator(object):

    typemap = {'Intron': 'intron'}
    fieldmap = {'sequence_id':  'seqid',
                'gene_type':    'type',
                'gene_product': 'product',
                'is_pseudo':    'pseudo'}

    def _assign_UTR_flanks(self, genes):
        logger.info('  assigning UTR directions (five_prime, three_prime)')
        for gene_id, gene in six.iteritems(genes):
            if not 'Gene Model' in gene:
                continue
            for transcript_id, genemodel in six.iteritems(gene['Gene Model']):
                logger.info('     processing transcript %s' % transcript_id)
                strand = None
                genemodel.sort(key=lambda feat: feat['Start'])
                for feature in genemodel:
                    if not strand:
                        strand = feature['Strand']
                    else:
                        if strand != feature['Strand']:
                            raise RuntimeError("inconsistent strands within " +
                                               "transcript %s" % transcript_id)
                seen_coding = False
                for feature in genemodel:
                    if feature['Type'] == 'UTR':
                        if seen_coding:
                            if strand == '+':
                                feature['Type'] = 'three_prime_UTR'
                            else:
                                feature['Type'] = 'five_prime_UTR'
                        else:
                            if strand == '+':
                                feature['Type'] = 'five_prime_UTR'
                            else:
                                feature['Type'] = 'three_prime_UTR'
                    elif feature['Type'] == 'CDS':
                        seen_coding = True

    def _process_gene_model(self, table, genes):
        #[Gene ID]       [Transcript ID(s)]      [Name]  [Type]  [Sequence_id]   [Start] [End]   [Is Reversed]
        for entry in table:
            entry['Gene ID'] = entry['Gene ID'].split(', ')[0]
            if not genes[entry['Gene ID']]:
                raise RuntimeError("orphan gene model with ID '%s' doesn't " +
                                   "have info in ByTaxon result" % entry['Gene ID'])
            target_gene = genes[entry['Gene ID']]
            target_transcripts = entry['Transcript ID(s)'].split(', ')
            assert(target_transcripts is not None)
            if entry['Type'] in self.typemap:
                entry['Type'] = self.typemap[entry['Type']]
            entry['Seqid'] = entry['Sequence_id']
            entry['Start'] = int(entry['Start'])
            entry['End'] = int(entry['End'])
            entry['Strand'] = "+-"[int(entry['Is Reversed'])]
            for tt in target_transcripts:
                if not 'Gene Model' in target_gene:
                    target_gene['Gene Model'] = {}
                if not tt in target_gene['Gene Model']:
                    target_gene['Gene Model'][tt] = []
                target_gene['Gene Model'][tt].append(entry)
        self._assign_UTR_flanks(genes)

    def _process_pubmed(self, table, genes):
        # [Gene ID]         [PubMed ID]     [doi]   [Title] [Authors]
        for entry in table:
            entry['Gene ID'] = entry['Gene ID'].split(', ')[0]
            if not genes[entry['Gene ID']]:
                raise RuntimeError("orphan gene model with ID '%s' doesn't " +
                                   "have info in ByTaxon result" % entry['Gene ID'])
            target_gene = genes[entry['Gene ID']]
            if not 'Dbxref' in target_gene:
                    target_gene['Dbxref'] = []
            target_gene['Dbxref'].append("PMID:%s" % str(entry['PubMed ID']))

    def _process_go(self, table, genes):
        # [Gene ID]         [Transcript ID(s)]      [Ontology]
        # [GO ID] [GO Term Name]  [Source]        [Evidence Code] [Is Not]
        # [Reference] [Evidence Code Support]
        for entry in table:
            entry['Gene ID'] = entry['Gene ID'].split(', ')[0]
            if not genes[entry['Gene ID']]:
                raise RuntimeError("orphan gene model with ID '%s' doesn't " +
                                   "have info in ByTaxon result" % entry['Gene ID'])
            target_gene = genes[entry['Gene ID']]
            target_transcripts = entry['Transcript ID(s)'].split(', ')
            assert(target_transcripts is not None)
            for tt in target_transcripts:
                if not 'GO Terms' in target_gene:
                    target_gene['GO Terms'] = {}
                if not tt in target_gene['GO Terms']:
                    target_gene['GO Terms'][tt] = []
                target_gene['GO Terms'][tt].append(entry)

    def _create_gene(self, genes, item):
        logger.info('  creating gene %s' % item['id'])
        v = {}
        genes[item['id']] = v
        v['ID'] = item['id']
        for f in item['fields']:
            if f['name'] in self.fieldmap:
                v[self.fieldmap[f['name']]] = f['value']
            elif f['name'] == 'gene_location_text':
                m = re.search('^([^:]+):([0-9,]+)\.\.([0-9,]+)\((.)\)', f['value'])
                if m:
                    v['start'] = int(m.group(2).replace(',', ''))
                    v['stop'] = int(m.group(3).replace(',', ''))
                    v['strand'] = m.group(4)
            else:
                v[f['name']] = f['value']

    def _get_json(self, url, params=None):
        logger.info('  retrieving %s' % url)
        s = get_session(self.baseurl, self.login)
        res = s.get(url, verify=True, params=params)
        if "autologin" in res.text or 'Login</title>' in res.text:
            raise RuntimeError("Login failed -- please check user credentials.")
        if(res.ok):
            return res.json()
        else:
            res.raise_for_status()

    def _parse_table(self, tablestring, tablename):
        logger.info('  parsing table %s' % tablename)
        s = StringIO(tablestring)
        i = 0
        headers = None
        out = []
        for line in s:
            if i == 0:
                headers = re.findall('\[([^[]+)\]', line)
                i = i + 1
            else:
                res = {}
                if not headers:
                    raise RuntimeError("table headers not found trying to " +
                                       "query table '%s'" % tablename)
                vals = line.split("\t")
                if len(vals) < len(headers):
                    raise RuntimeError("invalid number of columns (%s) " +
                                       "reading line '%s' in table '%s' " +
                                       "(expected %s)" % (str(len(vals)),
                                                          line, tablename,
                                                          str(len(headers))))
                j = 0
                for val in vals:
                    if j < len(headers):
                        res[headers[j]] = val
                    j = j + 1
                out.append(res)
        return out

    def _get_table(self, table):
        s = get_session(self.baseurl, self.login)

        query_payload = {
            "questionDefinition": {
                "questionName": "GeneQuestions.GenesByTaxonGene",
                "parameters": {"organism": self.organism},
                "viewFilters": [],
                "filters": []
            },
            "formatting": {
                "formatConfig": {
                    "tables": [table],
                    "includeHeader": True,
                    "attachmentType": "plain"
                },
                "format": "tableTabular"
            }
        }

        r = s.post(self.baseurl + "/service/answer",
                   data=json.dumps(query_payload),
                   headers={'Content-Type': 'application/json'})
        r.raise_for_status()
        return self._parse_table(r.text, table)

    def __init__(self, baseurl, organism, cache=False, login=None):
        self.baseurl = baseurl
        self.organism = organism
        self.login = parse_login(login)
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
                       'sequence_id']
        self.tables = ['GeneModelDump',
                       'Notes',
                       'PubMed',
                       'GOTerms']

        # get gene centric information
        params = {'organism': self.organism,
                  'o-fields': ','.join(self.fields)}

        genes = {}

        # get gene level data from JSON webservice
        url = '{0}/webservices/GeneQuestions/GenesByTaxon.json'.format(baseurl)
        taxon_json = self._get_json(url, params)
        for v in taxon_json['response']['recordset']['records']:
            self._create_gene(genes, v)

        # do fast table queries via new web service
        tbl = self._get_table('GeneModelDump')
        self._process_gene_model(tbl, genes)
        tbl = self._get_table('PubMed')
        self._process_pubmed(tbl, genes)
        tbl = self._get_table('GOTerms')
        self._process_go(tbl, genes)

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
                             'Exon': 'pseudogenic_exon',
                             'CDS': 'pseudogenic_exon'}

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
