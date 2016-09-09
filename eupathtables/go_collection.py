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

import warnings
import datetime
import six
try:
    from exceptions import RuntimeWarning
except:
    pass


class GOCollection(object):
    def __init__(self, taxon_id):
        self.gos = {}
        self.taxon_id = taxon_id
        self.aspects = {'Biological Process': 'P',
                        'Molecular Function': 'F',
                        'Cellular Component': 'C'}

    def _aspect2oneletter(self, aspect):
        return self.aspects[aspect]

    def add_item(self, v, source='EuPathDB'):
        """Add GO terms from an item parsed from EuPathTables to the set."""
        if 'GO Terms' in v:
            for go in v['GO Terms']:
                aspect = None
                object_symbol = v['ID']
                if 'name' in v:
                    object_symbol = v['name']
                try:
                    aspect = self._aspect2oneletter(go['Ontology'])
                except:
                    warnings.warn("couldn't get aspect for %s (%s)"
                                  % (go['GO ID'], go['GO Term Name']),
                                  RuntimeWarning)
                    continue
                if aspect and go['GO ID']:
                    self.add_generic(source, v['ID'], object_symbol, '',
                                     go["GO ID"], 'GO_REF:0000002',
                                     go['Evidence Code'], '', aspect,
                                     v['product'], '', 'gene', source)

    def add_from_gaf_iterator(self, it, stream=None):
        for item in it:
            # XXX this assumes object IDs are on protein level! Check if gene symbol is usable
            if (not stream) or (stream and (item['object_id'] in stream.uniprots)):
                if stream and (item['object_id'] in stream.uniprots):
                    item['object_id'] = stream.uniprots[item['object_id']]
                    item['object_type'] = 'gene'
                self.add_generic(item['db'], item['object_id'],
                                 item['object_symbol'], item['qualifier'],
                                 item['go_id'], item['dbref'],
                                 item['evidence_code'], item['withfrom'],
                                 item['aspect'], item['object_name'],
                                 item['object_synonym'], item['object_type'],
                                 item['assigned_by'])

    def add_generic(self, db='', object_id='', object_symbol='', qualifier='',
                    go_id='', dbref='', evidence_code='', withfrom='',
                    aspect='', object_name='', object_synonym='',
                    object_type='gene', assigned_by=''):
        """Add GO terms parsed separately to the set."""
        assert(go_id)
        new_item = {'db': db,  'object_id': object_id,
                    'object_symbol': object_symbol, 'qualifier': qualifier,
                    'go_id': go_id, 'dbref': dbref,
                    'evidence_code': evidence_code, 'withfrom': withfrom,
                    'aspect': aspect, 'object_name': object_name,
                    'object_synonym': object_synonym,
                    'object_type': object_type, 'taxon': self.taxon_id,
                    'assigned_by': assigned_by}
        if not go_id in self.gos:
            self.gos[go_id] = []
        self.gos[go_id].append(new_item)

    def to_gafv1(self, outstream):
        """Output the contents of this collection in GAF 1.0 format."""
        outstream.write("!gaf-version: 1.0\n")
        for k, v in six.iteritems(self.gos):
            for it in v:
                outstream.write(it['db'] + "\t" +
                                it['object_id'] + "\t" +
                                it['object_symbol'] + "\t" +
                                it['qualifier'] + "\t" +
                                it['go_id'] + "\t" +
                                it['dbref'] + "\t" +
                                it['evidence_code'] + "\t" +
                                it['withfrom'] + "\t" +
                                it['aspect'] + "\t" +
                                it['object_name'] + "\t" +
                                it['object_synonym'] + "\t" +
                                it['object_type'] + "\t" +
                                "taxon:" + str(it['taxon']) + "\t" +
                                datetime.date.today().strftime('%Y%m%d') + "\t" +
                                it['assigned_by'] + "\n")
