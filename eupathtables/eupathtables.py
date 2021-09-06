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

    typemap = {'Intron': 'intron',
               'Exon': 'exon'}
    fieldmap = {'sequence_id':  'seqid',
                'gene_type':    'type',
                'gene_product': 'product',
                'is_pseudo':    'pseudo'}

    def _assign_UTR_flanks(self, genes):
        logger.info('  assigning UTR directions (five_prime, three_prime)')
        for gene_id, gene in genes.items():
            if not 'Gene Model' in gene:
                continue
            for transcript_id, genemodel in gene['Gene Model'].items():
                logger.info('     processing transcript %s' % transcript_id)
                strands = set(f['Strand'] for f in genemodel)
                if len(strands) != 1:
                    raise RuntimeError("inconsistent strands within " +
                                       "transcript %s" % transcript_id)
                strand = strands.pop()
                genemodel.sort(key=lambda feat: feat['Start'])                            
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
            entry['Gene ID'] = entry['Gene ID'].split(',')[0].strip()
            if entry['Gene ID'] not in genes:
                raise RuntimeError("orphan gene model with ID '%s' doesn't " +
                                "have info in ByTaxon result" % entry['Gene ID'])
            target_gene = genes[entry['Gene ID']]
            if not 'Gene Model' in target_gene:
                target_gene['Gene Model'] = {}
            if entry['Type'] in self.typemap:
                entry['Type'] = self.typemap[entry['Type']]
            entry['Seqid'] = entry['Sequence_id']
            entry['Start'] = int(entry['Start'])
            entry['End'] = int(entry['End'])
            entry['Strand'] = "+-"[int(entry['Is Reversed'])]
            for transcript in entry['Transcript ID(s)'].split(','):
                target_transcript = transcript.strip()
                assert target_transcript
                if not target_transcript in target_gene['Gene Model']:
                    target_gene['Gene Model'][target_transcript] = []
                target_gene['Gene Model'][target_transcript].append(entry)
        self._assign_UTR_flanks(genes)

    def _process_pubmed(self, table, genes):
        # [Gene ID]         [PubMed ID]     [doi]   [Title] [Authors]
        for entry in table:
            entry['Gene ID'] = entry['Gene ID'].split(',')[0].strip()
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
            entry['Gene ID'] = entry['Gene ID'].split(',')[0].strip()
            if entry['Gene ID'] not in genes:
                    raise RuntimeError("orphan gene model with ID '%s' doesn't " +
                                       "have info in ByTaxon result"%(entry['Gene ID']))
            target_gene = genes[entry['Gene ID']]
            if 'GO Terms' not in target_gene:
                target_gene['GO Terms'] = {}
            for transcript in entry['Transcript ID(s)'].split(','):
                target_transcript = transcript.strip()
                assert target_transcript                    
                if target_transcript not in target_gene['GO Terms']:
                    target_gene['GO Terms'][target_transcript] = []
                target_gene['GO Terms'][target_transcript].append(entry)
    
    def _process_gene_transcripts(self, table, genes):
        # [Gene ID] [Transcript]    [# exons]	
        # [Transcript length]   [Protein length]    [Transcript Type]
        for entry in table:            
            entry['Gene ID'] = entry['Gene ID'].split(',')[0].strip()
            if entry['Gene ID'] not in genes:
                raise RuntimeError("orphan gene model with ID '%s' doesn't " +
                                    "have info in ByTaxon result"%(entry['Gene ID']))
            target_gene = genes[entry['Gene ID']]            
            target_gene['transcript_type'] = entry['Transcript Type']
            target_gene['transcript_length'] = entry['Transcript length']

    def _create_gene(self, genes, item):
        v = {}
        geneId = item['displayName']
        logger.info('   creating gene %s' % geneId)
        genes[geneId] = v
        v['ID'] = geneId
        for key, val in item['attributes'].items():
            if key in self.fieldmap:
                v[self.fieldmap[key]] = val
            elif key == 'gene_location_text':
                m = re.search('^([^:]+):([0-9,]+)\.\.([0-9,]+)\((.)\)', val)
                if m:
                    v['start'] = int(m.group(2).replace(',', ''))
                    v['stop'] = int(m.group(3).replace(',', ''))
                    v['strand'] = m.group(4)
            else:
                v[key] = val

    def _get_json(self):
        url = self.url + '/standard'
        logger.info('  retrieving %s' % url)        

        params = {
            'organism': self.organism.replace('#',"%23"),
            'reportConfig': {
                "attributes": self.fields,
            }
        }

        # workaround for problematic urlencoding if organism contains '#'
        params = "?organism=" + params['organism'] + \
                 "&reportConfig=" + json.dumps(params['reportConfig'])
        res = self.session.get(url + params, verify=True)
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
                headers = line.split('\t')
                i += 1
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
                    j += 1
                out.append(res)
        return out

    def _get_table(self, table):
        url = self.url + '/tableTabular'
        logger.info(self.organism)

        query_payload = {
            "searchConfig": {
                "parameters": {"organism": self.organism},
            },
            "reportConfig": {
                "tables": [table],
                "includeHeader": True,
                "attachmentType": "plain",
            }
        }

        r = self.session.post(url, data=json.dumps(query_payload), headers={'Content-Type': 'application/json'})
        r.raise_for_status()
        return self._parse_table(r.text, table)

    def __init__(self, baseurl, organism, cache=False, session=None):
        self.baseurl = baseurl
        self.organism = organism
        self.session = session
        self.fields = ['annotated_go_function',
                       'gene_location_text',
                       'location_text',
                       'gene_type',
                       'gene_product',
                       'primary_key',
                       'ec_numbers',
                       'chromosome',
                       'gene_type',
                       'is_pseudo',
                       'gene_name',
                       'uniprot_id',
                       'sequence_id']

        self.url = '{0}/service/record-types/transcript/searches/GenesByTaxon/reports'.format(baseurl)

        genes = {}

        # get gene centric information
        taxon_json = self._get_json()
        for v in taxon_json['records']:
            self._create_gene(genes, v)

        # do fast table queries via new web service
        tbl = self._get_table('GeneModelDump')
        self._process_gene_model(tbl, genes)
        tbl = self._get_table('PubMed')
        self._process_pubmed(tbl, genes)
        tbl = self._get_table('GOTerms')
        self._process_go(tbl, genes)
        tbl = self._get_table('GeneTranscripts')
        self._process_gene_transcripts(tbl, genes)

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
            self.pseudo_typemaps = {'gene': 'pseudogene',
                                    'mRNA': 'pseudogenic_transcript',
                                    'rRNA': 'pseudogenic_transcript',
                                    'tRNA': 'pseudogenic_transcript',
                                    'exon': 'pseudogenic_exon',
                                    'CDS': 'pseudogenic_exon'}

        def finaltype(self, v, mtype):
            if 'pseudo' in v and v['pseudo'] == 'Yes':
                if mtype in self.pseudo_typemaps:
                    return self.pseudo_typemaps[mtype]
            return mtype

        def make_noncoding(self, v, gene):
            if 'Gene Model' in v:
                for transcript_id, transcript_model in six.iteritems(v['Gene Model']):                    
                    transcript = gt.extended.FeatureNode.create_new(v['seqid'], self.finaltype(v, v['transcript_type']),
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
                    if 'gene_name' in v and v['gene_name'] is not None:
                        gene.add_attribute("Name", v['gene_name'])

                    # track UniProtID -> gene mappings
                    if 'uniprot_id' in v and v['uniprot_id'] is not None:
                        for u in v['uniprot_id'].split(','):
                            if len(u) > 0:
                                self.uniprots[u] = v['ID']                    

                    # protein coding gene
                    if v['type'] == 'protein coding gene':
                        self.make_coding(v, gene)
                    # non-coding RNA
                    elif v['type'] == 'ncRNA gene':
                        self.make_noncoding(v, gene)
                    
                    break

                #except:
                #    if v:
                #        sys.stderr.write("error creating feature for %s\n" % v['ID'])
                #    else:
                #        sys.stderr.write("error creating feature , no ID yet\n")
                #    continue

            return gene
