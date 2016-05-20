#!/usr/bin/env python
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

import exceptions

class GAFError(exceptions.RuntimeError):
    pass

class GAFIterator(object):
    def __init__(self, instream):
        self.instream = instream
        self.keys = ['db', 'object_id', 'object_symbol', 'qualifier',
                     'go_id', 'dbref', 'evidence_code', 'withfrom',
                     'aspect', 'object_name', 'object_synonym',
                     'object_type', 'taxon', 'datetime', 'assigned_by']

    def __iter__(self):
        return self

    def next(self):
        return self.__next__()

    def __next__(self):
        l = None
        while not l or len(l) == 0 or l[0] == '!':
            l = next(self.instream)
        if l and len(l) > 0:
            arr = l.split("\t")
            if len(arr) < 15:
                raise GAFError("less than 15 columns in line: %s" % (l))
            out = dict(zip(self.keys, arr))
            return out
        else:
            raise exceptions.StopIteration
