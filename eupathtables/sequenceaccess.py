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
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class SequenceProvider(object):

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

    def __init__(self, baseurl, organism, login=None):
        self.login = parse_login(login)
        self.baseurl = baseurl
        url = ('{0}/webservices/GenomicSequenceQuestions/' +
               'SequencesByTaxon.json' +
               '?organism={1}&o-fields=primary_key').format(baseurl, organism.replace('#', "%23"))
        res = self._get_json(url)
        seqids = []
        for v in res['response']['recordset']['records']:
            seqids.append(v['id'])
        logger.info('  needing to retrieve %s seqs for organism %s' % (len(seqids), organism))
        parsed_url = urlparse(baseurl)
        self.baseurl = "%s://%s" % (parsed_url.scheme, parsed_url.netloc)
        # XXX: hard coded FungiDB, please make configurable!
        payload = {'project_id': 'FungiDB',
                   'ids': '\n'.join(seqids)}
        url = ('{0}/cgi-bin/contigSrt').format(self.baseurl)
        s = get_session(self.baseurl, self.login)
        res = s.post(url, data=payload)
        if not res.ok:
            res.raise_for_status()
        self.out = res.text

    def to_file(self, filename):
        f = open(filename, 'w+')
        f.write(self.out)
        f.close()
