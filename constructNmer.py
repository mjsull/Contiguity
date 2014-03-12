import sys
import array

class nmerDict:
    def __init__(self, readfile, nmersize):
        self.readfile = readfile
        self.nmersize = nmersize
    def get_nmer_freq(self):
        nmersize, reads  = self.nmersize, self.readfile
        nucl = set('atcg')
        self.nmerdict = {}
        if reads[-3:] == '.gz':
            readfile = gzip.open(reads)
        else:
            readfile = open(reads)
        getseq = False
        getfaseq = False
        for line in readfile:
            if line.startswith('@'):
                getseq = True
            elif line.startswith('>'):
                if getfaseq:
                    for i in range(0, len(seq) - nmersize):
                        nmer = seq[i:i+nmersize]
                        if nmer[nmersize/2] == 'a' or nmer[nmersize/2] == 'c':
                            nmer = nmer[::-1]
                            nmer = nmer.translate(transtab)
                        if nmer in self.nmerdict:
                            self.nmerdict[nmer] += 1
                        else:
                            if set(nmer).issubset(nucl):
                                self.nmerdict[nmer] = 1
                getfaseq = True
                seq = ''
            elif getfaseq:
                seq += line.rstrip().lower()
            elif getseq:
                getseq = False
                seq = line.rstrip().lower()
                for i in range(0, len(seq) - nmersize + 1):
                    nmer = seq[i:i+nmersize]
                    if nmer[nmersize/2] == 'a' or nmer[nmersize/2] == 'c':
                        nmer = nmer[::-1]
                        nmer = nmer.translate(transtab)
                    if nmer in self.nmerdict:
                        self.nmerdict[nmer] += 1
                    else:
                        if set(nmer).issubset(nucl):
                            self.nmerdict[nmer] = 1
        if getfaseq:
            for i in range(0, len(seq) - nmersize + 1):
                nmer = seq[i:i+nmersize]
                if nmer[nmersize/2] == 'a' or nmer[nmersize/2] == 'c':
                    nmer = nmer[::-1]
                    nmer = nmer.translate(transtab)
                if nmer in self.nmerdict:
                    self.nmerdict[nmer] += 1
                else:
                    if set(nmer).issubset(nucl):
                        self.nmerdict[nmer] = 1
        out = open(self.workingDir.get() + '/nmers', 'w')
        for i in self.nmerdict:
            out.write(i + '\t' + str(self.nmerdict[i]) + '\n')
        out.close()
        out = open(self.workingDir.get() + '/nmerfreq.csv', 'w')
        freqDict = {}
        for i in self.nmerdict:
            if self.nmerdict[i] in freqDict:
                freqDict[self.nmerdict[i]] += 1
            else:
                freqDict[self.nmerdict[i]] = 1
        for i in range(1, max(freqDict) + 1):
            if i in freqDict:
                out.write(str(i) + '\t' + str(freqDict[i]) + '\n')
            else:
                out.write(str(i) + '\t0\n')
        out.close()
        self.nmerfile = self.workingDir.get() + '/nmers'
        return False # thread has not been aborted

    def get_nmer_freq(self):
        nmersize, reads  = self.nmersize.get(), self.readfile.get()
        nucl = set('atcg')
        self.nmerdict = {}
        if reads[-3:] == '.gz':
            readfile = gzip.open(reads)
        else:
            readfile = open(reads)
        getseq = False
        getfaseq = False
        for line in readfile:
            if not self.abortqueue.empty():
                if self.aborttime():
                    return True
            if line.startswith('@'):
                getseq = True
            elif line.startswith('>'):
                if getfaseq:
                    for i in range(0, len(seq) - nmersize):
                        nmer = seq[i:i+nmersize]
                        if nmer[nmersize/2] == 'a' or nmer[nmersize/2] == 'c':
                            nmer = nmer[::-1]
                            nmer = nmer.translate(transtab)
                        if nmer in self.nmerdict:
                            self.nmerdict[nmer] += 1
                        else:
                            if set(nmer).issubset(nucl):
                                self.nmerdict[nmer] = 1
                getfaseq = True
                seq = ''
            elif getfaseq:
                seq += line.rstrip().lower()
            elif getseq:
                getseq = False
                seq = line.rstrip().lower()
                for i in range(0, len(seq) - nmersize + 1):
                    nmer = seq[i:i+nmersize]
                    if nmer[nmersize/2] == 'a' or nmer[nmersize/2] == 'c':
                        nmer = nmer[::-1]
                        nmer = nmer.translate(transtab)
                    if nmer in self.nmerdict:
                        self.nmerdict[nmer] += 1
                    else:
                        if set(nmer).issubset(nucl):
                            self.nmerdict[nmer] = 1
        if getfaseq:
            for i in range(0, len(seq) - nmersize + 1):
                nmer = seq[i:i+nmersize]
                if nmer[nmersize/2] == 'a' or nmer[nmersize/2] == 'c':
                    nmer = nmer[::-1]
                    nmer = nmer.translate(transtab)
                if nmer in self.nmerdict:
                    self.nmerdict[nmer] += 1
                else:
                    if set(nmer).issubset(nucl):
                        self.nmerdict[nmer] = 1
        out = open(self.workingDir.get() + '/nmers', 'w')
        for i in self.nmerdict:
            out.write(i + '\t' + str(self.nmerdict[i]) + '\n')
        out.close()
        out = open(self.workingDir.get() + '/nmerfreq.csv', 'w')
        freqDict = {}
        for i in self.nmerdict:
            if self.nmerdict[i] in freqDict:
                freqDict[self.nmerdict[i]] += 1
            else:
                freqDict[self.nmerdict[i]] = 1
        for i in range(1, max(freqDict) + 1):
            if i in freqDict:
                out.write(str(i) + '\t' + str(freqDict[i]) + '\n')
            else:
                out.write(str(i) + '\t0\n')
        out.close()
        self.nmerfile = self.workingDir.get() + '/nmers'
        return False # thread has not been aborted