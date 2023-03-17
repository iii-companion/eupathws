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

import os
import sys
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)
from eupathtables.eupathtables import WebServiceIterator, TableInStream
from eupathtables.sequenceaccess import SequenceProvider
from gt.extended import GFF3OutStream
from contextlib import contextmanager
import ctypes
import io
import json
import tempfile
import requests
libc = ctypes.CDLL(None)
c_stdout = ctypes.c_void_p.in_dll(libc, 'stdout')


@contextmanager
def stdout_redirector(stream):
    # The original fd stdout points to. Usually 1 on POSIX systems.
    original_stdout_fd = sys.stdout.fileno()

    def _redirect_stdout(to_fd):
        """Redirect stdout to the given file descriptor."""
        # Flush the C-level buffer stdout
        libc.fflush(c_stdout)
        # Flush and close sys.stdout - also closes the file descriptor (fd)
        sys.stdout.close()
        # Make original_stdout_fd point to the same file as to_fd
        os.dup2(to_fd, original_stdout_fd)
        # Create a new sys.stdout that points to the redirected fd
        sys.stdout = io.TextIOWrapper(os.fdopen(original_stdout_fd, 'wb'))

    # Save a copy of the original stdout fd in saved_stdout_fd
    saved_stdout_fd = os.dup(original_stdout_fd)
    try:
        # Create a temporary file and redirect stdout to it
        tfile = tempfile.TemporaryFile(mode='w+b')
        _redirect_stdout(tfile.fileno())
        # Yield to caller, then redirect stdout back to the saved fd
        yield
        _redirect_stdout(saved_stdout_fd)
        # Copy contents of temporary file to the given stream
        tfile.flush()
        tfile.seek(0, io.SEEK_SET)
        stream.write(tfile.read())
    finally:
        tfile.close()
        os.close(saved_stdout_fd)


def has_all_files(organism):
    for sfx in ['_Genome.fasta', '.gaf', '.gff3', '_Proteins.fasta']:
        if not os.path.isfile("%s%s" % (organism, sfx)):
            return False
    return True


def download_file(url, local_filename=None):
    if not local_filename:
        local_filename = url.split('/')[-1]
    r = requests.get(url, stream=True)
    with open(local_filename, 'wb') as f:
        for chunk in r.iter_content(chunk_size=1024):
            if chunk:
                f.write(chunk)
    return local_filename


def create_fileset_for_organism(baseurl, organism, session):
    logger.info("fetching data for organism %s" % organism)

    wsit = WebServiceIterator(baseurl, organism, cache=False,
                              session=session)
    laststream = table_in_stream = TableInStream(wsit)
    laststream = GFF3OutStream(table_in_stream)

    # capture stdout from output stream
    f = io.BytesIO()
    with stdout_redirector(f):
        fn = laststream.next_tree()
        while fn:
            fn = laststream.next_tree()

    # write out GFF for organism
    with open('%s.gff3' % organism.replace('/', '_'), 'w') as fd:
        fd.write("{0}".format(f.getvalue().decode('utf-8')))

    # write out chromosomes JSON for organism
    with open('%s_chr.json' % organism.replace('/', '_'), 'w') as fp:
        json.dump(wsit.chromosomes, fp)

    # write out metadata JSON for organism
    with open('%s_meta.json' % organism.replace('/', '_'), 'w') as fp:
        json.dump(wsit.meta, fp)

    # download FASTA for organism
    sp = SequenceProvider(baseurl, organism, session=session)
    sp.to_file("%s_Genome.fasta" % organism.replace('/', '_'))

    # download Proteins FASTA for organism
    sp = SequenceProvider(baseurl, organism, session=session, typ="protein")
    sp.to_file("%s_Proteins.fasta" % organism.replace('/', '_'))

    # write out GAF for organism
    table_in_stream.go_coll.to_gafv1(open('%s.gaf' % organism.replace('/', '_'), "w+"))
