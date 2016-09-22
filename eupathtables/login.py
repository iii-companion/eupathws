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

import requests
import logging
logger = logging.getLogger(__name__)


def get_session(baseurl, login=None):
    logger.info('  creating session')
    s = requests.session()
    if login:
        logger.info('  logging in')
        s.get(baseurl)
        login_payload = {
            'username': login.username,
            'password': login.password,
            "submit": "Login"
        }
        s.post("https://eupathdb.org/auth/bin/login", data=login_payload)
        s.get(baseurl)
    return s

def parse_login(login):
    if not login:
        return None
    vals = str(login).split(':')
    if len(vals) != 2:
        raise RuntimeError("login string '%s' does not have " +
                           "expected format user:password" % login)
    else:
        from collections import namedtuple
        Login = namedtuple("Login", "username password")
        return Login(vals[0], vals[1])
