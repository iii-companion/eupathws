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
import json
from eupathtables.login import parse_login, get_session


class GeneMetrics(object):

    def __init__(self, baseurl, session):
        self.baseurl = baseurl
        self.fields = ["organism",
                       "ncbi_tax_id",
                       "is_reference_strain",
                       "is_annotated_genome",
                       "URLGenomeFasta",
                       "URLgff"]
        url = "{0}/service/record-types/organism/searches/GeneMetrics/reports/standard".format(baseurl)
        
        query_payload = {
            "searchConfig": {
                "parameters": {},
            },
            "reportConfig": {
                "tables": [],
                "attributes": self.fields,
            }
        }

        res = session.post(url, data=json.dumps(query_payload), headers={'Content-Type': 'application/json'})
        self.orgs = collections.deque()
        if(res.ok):
            j = res.json()
        else:
            res.raise_for_status()
        for rec in j['records']:
            rec['attributes'].update({'displayName': rec['displayName']})
            self.orgs.append(rec['attributes'])

    def __iter__(self):
        return self

    def next(self):
        return self.__next__()

    def __next__(self):
        if len(self.orgs) > 0:
            return self.orgs.popleft()
        else:
            raise StopIteration
