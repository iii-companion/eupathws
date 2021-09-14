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

from bs4 import BeautifulSoup
from urllib.parse import urlparse
import logging
import urllib

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class SequenceProvider(object):

    def _get_json(self):
        logger.info('  retrieving %s' % self.url)
        
        res = self.session.get(self.url, verify=True)
        if "autologin" in res.text or 'Login</title>' in res.text:
            raise RuntimeError("Login failed -- please check user credentials.")
        if res.ok:
            return res.json()
        else:
            res.raise_for_status()

    def __init__(self, baseurl, organism, session=None):
        self.session = session
        field = "URLGenomeFasta"
        parsed_baseurl = urlparse(baseurl)
        baseurl = "%s://%s" % (parsed_baseurl.scheme, parsed_baseurl.netloc)
        self.url = '{0}/a/service/record-types/organism/searches/GenomeDataTypes/reports/standard?reportConfig={{\"attributes\":[\"{1}\"]}}'.format(baseurl, field)

        res = self._get_json()                

        self.download_url = ''
        for v in res['records']:
            # Remove html formatting from displayName
            if BeautifulSoup(v['displayName'], 'html.parser').text == organism:
                self.download_url = v['attributes'][field]
                break
        if not self.download_url:
            raise RuntimeError("Failed to establish fasta download location for organism %s" % organism)

    def to_file(self, filename):
        logger.info("Retrieving fasta file from %s" % self.download_url)
        try:
            urllib.request.urlretrieve(self.download_url, filename)
        except urllib.error.URLError as e:
            logging.error("Cannot retrieve file from url: {0}. Please check the URL is correct. In case of an outage at EuPathDB please try again later.\n\n{1}".format(self.download_url, e))
            raise SystemExit()
