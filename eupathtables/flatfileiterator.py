class FlatFileIterator(object):

    def __init__(self, instream):
        self.instream = instream
        self.eof = False
        self.gene = dict()

    def __iter__(self):
        return self

    def next(self):
        return self.__next__()

    def __next__(self):
        self.gene = dict()
        while True:
            l = self.instream.readline()
            if l == '' or l is None:
                if not self.eof:
                    self.eof = True
                    return self.gene
                else:
                    raise exceptions.StopIteration
            m = re.match(r"Gene ID: (.*)$", l)
            if m and m.group(1):
                self.gene['ID'] = m.group(1)
                continue
            m = re.match(r"Gene Type: (.*)$", l)
            if m and m.group(1):
                self.gene['type'] = m.group(1)
                continue
            m = re.match(r"UniProt ID: (.*)$", l)
            if m and m.group(1) and m.group(1) != 'null':
                self.gene['uniprot_id'] = m.group(1)
                continue
            m = re.match(r"Gene Name: (.*)$", l)
            if m and m.group(1) and m.group(1) != 'null':
                self.gene['name'] = m.group(1)
                continue
            m = re.match(r"Is Pseudo: (.*)$", l)
            if m and m.group(1):
                if m.group(1) == 'Yes':
                    self.gene['pseudo'] = True
                else:
                    self.gene['pseudo'] = False
                continue
            m = re.match(r"Annotated ([35]). UTR length: ([0-9]+)", l)
            if m and m.group(1) and m.group(2):
                self.gene["utr_%s" % m.group(1)] = int(m.group(2))
                continue
            m = re.match(r"Product Description: (.*)$", l)
            if m and m.group(1):
                self.gene['product'] = m.group(1)
                continue
            m = re.match(r"Genomic Location: (?P<seqid>[^ ]+): " +
                         "(?P<start>[0-9,]+) - (?P<stop>[0-9,]+) " +
                         "\((?P<strand>.)\)", l)
            if m:
                for k, v in m.groupdict().iteritems():
                    self.gene[k] = v.replace(',', '')
                continue
            m = re.match(r"TABLE: (.*)$", l)
            if m and m.group(1):
                thistable = []
                tablename = m.group(1)
                l = self.instream.readline()
                listidx = re.findall("\[([^\]]+)\]", l)
                l = self.instream.readline()
                while l != '\n' and l != '':
                    i = 0
                    thisline = dict()
                    for v in re.split('\t', l):
                        thisline[listidx[i]] = v.rstrip()
                        i = i + 1
                    thistable.append(thisline)
                    l = self.instream.readline()
                self.gene[tablename] = thistable
                continue
            if re.match('------------------', l):
                return self.gene