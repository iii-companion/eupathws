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

from eupathtables.login import parse_login, get_session
from urllib.parse import urlparse
import json
import logging
import re
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class SequenceProvider(object):

    def _get_json(self):
        url = self.url + '/standard'
        logger.info('  retrieving %s' % url)

        params = {
            "organism": self.organism.replace("#","%23"),
            "reportConfig": {
                "attributes": [
                    "primary_key",
                    "organism",
                    "formatted_length"
                ],
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

    def __init__(self, baseurl, organism, session=None):
        self.baseurl = baseurl
        self.organism = organism
        self.session = session
        
        self.url = '{0}/service/record-types/genomic-sequence/searches/SequencesByTaxon/reports'.format(baseurl)
        
        res = self._get_json()
        seqids = []
        for v in res['records']:
            seqids.append(v['displayName'])
        logger.info('  needing to retrieve %s seqs for organism %s' % (len(seqids), organism))
        parsed_url = urlparse(baseurl)
        self.baseurl = "%s://%s" % (parsed_url.scheme, parsed_url.netloc)
        project_id = parsed_url.netloc.split('.')[0].capitalize()
        project_id = re.sub('db$', 'DB', project_id)
        project_id = project_id.replace('Tritryp', 'TriTryp')
        project_id = project_id.replace('Vectorbase', 'VectorBase')
        payload = {'project_id': project_id,
                   'ids': '\n'.join(seqids)}
        url = ('{0}/cgi-bin/contigSrt').format(self.baseurl)
        res = self.session.post(url, data=payload)
        if not res.ok:
            res.raise_for_status()
        self.out = res.text

    def to_file(self, filename):
        f = open(filename, 'w+')
        f.write(self.out)
        f.close()
