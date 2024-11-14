from KarkSand import direct_kark_sort
import time
import pickle


class FmCheckpoints(object):
    ''' Manages rank checkpoints and handles rank queries, which are
        O(1) time, with the checkpoints taking O(m) space, where m is
        length of text. '''

    def __init__(self, bw, cpIval=1):
        ''' Scan BWT, creating periodic checkpoints as we go '''
        self.cps = {}        # checkpoints
        self.cpIval = cpIval # spacing between checkpoints
        tally = {}           # tally so far
        # Create an entry in tally dictionary and checkpoint map for
        # each distinct character in text
        for c in bw:
            if c not in tally:
                tally[c] = 0
                self.cps[c] = []
        # Now build the checkpoints
        for i, c in enumerate(bw):
            tally[c] += 1 # up to *and including*
            if i % cpIval == 0:
                for c in tally.keys():
                    self.cps[c].append(tally[c])

    def rank(self, bw, c, row):
        ''' Return # c's there are in bw up to and including row '''
        if row < 0 or c not in self.cps:
            return 0
        i, nocc = row, 0
        # Always walk to left (up) when calculating rank
        while (i % self.cpIval) != 0:
            if bw[i] == c:
                nocc += 1
            i -= 1
        return self.cps[c][i // self.cpIval] + nocc


class FmIndex():
    ''' O(m) size FM Index, where checkpoints and suffix array samples are
        spaced O(1) elements apart.  Queries like count() and range() are
        O(n) where n is the length of the query.  Finding all k
        occurrences of a length-n query string takes O(n + k) time.

        Note: The spacings in the suffix array sample and checkpoints can
        be chosen differently to achieve different bounds. '''

    @staticmethod
    def downsampleSuffixArray(sa, n=4):
        ''' Take only the suffix-array entries for every nth suffix.  Keep
            suffixes at offsets 0, n, 2n, etc with respect to the text.
            Return map from the rows to their suffix-array values. '''
        ssa = {}
        for i, suf in enumerate(sa):
            # We could use i % n instead of sa[i] % n, but we lose the
            # constant-time guarantee for resolutions
            if suf % n == 0:
                ssa[i] = suf
        return ssa

    def __init__(self, t, cpIval=1, ssaIval=4):
        if t[-1] != '$':
            t += '$' # add dollar if not there already
        # Get BWT string and offset of $ within it
        start_time = time.time()
        sa = direct_kark_sort(t)
        step_time = (time.time() - start_time)/60
        print(f'SA time: {step_time} min')

        self.bwt, self.dollarRow = self.bwtFromSa(t, sa)
        # Get downsampled suffix array, taking every 1 out of 'ssaIval'
        # elements w/r/t T
        self.ssa = self.downsampleSuffixArray(sa, ssaIval)
        self.slen = len(self.bwt)
        # Make rank checkpoints
        self.cps = FmCheckpoints(self.bwt, cpIval)
        # Calculate # occurrences of each character
        tots = dict()
        for c in self.bwt:
            tots[c] = tots.get(c, 0) + 1
        # Calculate concise representation of first column
        self.first = {}
        totc = 0
        for c, count in sorted(tots.items()):
            self.first[c] = totc
            totc += count

    def bwtFromSa(self, t, sa=None):
        ''' Given T, returns BWT(T) by way of the suffix array. '''
        bw = []
        dollarRow = None
        if sa is None:
            sa = direct_kark_sort(t)
        for si in sa:
            if si == 0:
                dollarRow = len(bw)
                bw.append('$')
            else:
                bw.append(t[si-1])
        return ''.join(bw), dollarRow # return string-ized version of list bw

    def count(self, c):
        ''' Return number of occurrences of characters < c '''
        if c not in self.first:
            # (Unusual) case where c does not occur in text
            for cc in sorted(self.first.iterkeys()):
                if c < cc: return self.first[cc]
            return self.first[cc]
        else:
            return self.first[c]

    def range(self, p):
        ''' Return range of BWM rows having p as a prefix '''
        l, r = 0, self.slen - 1 # closed (inclusive) interval
        for i in range(len(p)-1, -1, -1): # from right to left
            l = self.cps.rank(self.bwt, p[i], l-1) + self.count(p[i])
            r = self.cps.rank(self.bwt, p[i], r)   + self.count(p[i]) - 1
            if r < l:
                break
        return l, r+1

    def resolve(self, row):
        ''' Given BWM row, return the original position of a character in this row
        with respect to the original input text T.'''
        def stepLeft(row):
            ''' Step left according to character in given BWT row '''
            c = self.bwt[row]
            return self.cps.rank(self.bwt, c, row-1) + self.count(c)
        nsteps = 0
        while row not in self.ssa:
            row = stepLeft(row)
            nsteps += 1
        return self.ssa[row] + nsteps

    def hasSubstring(self, p):
        ''' Return true if and only if p is substring of indexed text '''
        l, r = self.range(p)
        return r > l

    def hasSuffix(self, p):
        ''' Return true if and only if p is suffix of indexed text '''
        l, r = self.range(p)
        off = self.resolve(l)
        return r > l and off + len(p) == self.slen-1

    def occurrences(self, p):
        ''' Return offsets for all occurrences of p, in no particular order '''
        l, r = self.range(p)
        return [ self.resolve(x) for x in range(l, r) ]


    def hasSubstringSpecial(self, p):
        l, r = 0, self.slen - 1  # closed interval
        intervals = []
        for i in range(len(p)-1, -1, -1):  # from right to left
            if p[i] != '#':
                if intervals:
                    new_intervals = []
                    for l, r in intervals:
                        l = self.cps.rank(self.bwt, p[i], l-1) + self.count(p[i])
                        r = self.cps.rank(self.bwt, p[i], r) + self.count(p[i]) - 1
                        if l <= r:
                            new_intervals.append((l, r))
                    intervals = new_intervals

                else:
                    l = self.cps.rank(self.bwt, p[i], l-1) + self.count(p[i])
                    r = self.cps.rank(self.bwt, p[i], r) + self.count(p[i]) - 1

            else:
                visited = set()  # process each character once
                for j in range(l, r+1): # calculate new left and right boundaries for each character in current range
                    c = self.bwt[j]
                    if c not in visited:
                        visited.add(c)
                        l_c = self.cps.rank(self.bwt, c, l-1) + self.count(c)
                        r_c = self.cps.rank(self.bwt, c, r) + self.count(c) - 1
                        if l_c <= r_c:  # valid interval
                            intervals.append((l_c, r_c))

        for l, r in intervals:
            if l <= r:
                return True
        return False
    
    def save_to_file(self, filename):
        with open(filename, 'wb') as file:
            pickle.dump(self, file)
        print(f'FmIndex instance saved to {filename}')

