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
import requests


class GeneMetrics(object):

    def __init__(self, baseurl, authtkt=None):
        self.baseurl = baseurl
        self.fields = ["organism",
                       "ncbi_tax_id",
                       "is_reference_strain",
                       "is_annotated_genome",
                       "URLGenomeFasta",
                       "URLgff"]
        url = ('{0}/webservices/OrganismQuestions/GeneMetrics.json?o-fields={1}').format(self.baseurl, ','.join(self.fields))
        if authtkt:
            url = url + ("&auth_tkt=%s" % (str(authtkt)))
        print(url)
        res = requests.get(url, verify=True)
        print(res.text)
        self.orgs = collections.deque()
        if(res.ok):
            j = res.json()
        else:
            res.raise_for_status()
        for v in j['response']['recordset']['records']:
            item = {}
            for f in v['fields']:
                item[f['name']] = f['value']
            self.orgs.append(item)

    def __iter__(self):
        return self

    def next(self):
        return self.__next__()

    def __next__(self):
        if len(self.orgs) > 0:
            return self.orgs.popleft()
        else:
            raise StopIteration
