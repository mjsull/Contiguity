# Contiguity - Tool for visualising assemblies
# Copyright (C) 2013-2015 Mitchell Sullivan
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Mitchell Sullivan
# mjsull@gmail.com
# School of Chemistry & Molecular Biosciences
# The University of Queensland
# Brisbane, QLD 4072.
# Australia


import sys, subprocess, string

transtab = string.maketrans('atcgATCG', 'tagcTAGC')
theident = 85


class contig:
    def __init__(self, name, forseq, revseq):
        self.name = name
        self.forseq = forseq
        self.revseq = revseq
        self.to = []
        self.fr = []
        self.length = len(forseq)

contigfile = open(sys.argv[1])
getline = 0
contigDict = {}
edgelist = []
templg = open('templg.fa', 'w')
bignodes = set()
count = 0
first = True
temped = open('tempedge.fa', 'w')
for line in contigfile:
    if line.startswith('NODE'):
        name = line.split()[1]
        seq = line.split()[-1]
        revseq = seq[::-1].translate(transtab)
        aninstance = contig(name, seq, revseq)
        contigDict[name] = aninstance
        templg.write('>' + name + '\n' + seq + '\n')
    elif line.startswith('EDGE'):
        junk, a, dira, b, dirb, overlap = line.split()
        if overlap == '.':
            overlap = ''
        if not 'n' in overlap:
            temped.write('>' + str(count) + '\n')
            if dira == 'True':
                if dirb == 'True':
                    if overlap.isdigit():
                        contigDict[a].to.append((b, True, int(overlap)))
                        contigDict[b].fr.append((a, False, int(overlap)))
                        temped.write(contigDict[a].forseq + contigDict[b].forseq[int(overlap):] + '\n')
                        edgelist.append((a, b, len(contigDict[a].forseq + contigDict[b].forseq[int(overlap):])))
                    else:
                        contigDict[a].to.append((b, True, overlap))
                        contigDict[b].fr.append((a, False, overlap[::-1].translate(transtab)))
                        temped.write(contigDict[a].forseq + overlap + contigDict[b].forseq + '\n')
                        edgelist.append((a, b, len(contigDict[a].forseq + overlap + contigDict[b].forseq)))
                else:
                    if overlap.isdigit():
                        contigDict[a].to.append((b, False, int(overlap)))
                        contigDict[b].to.append((a, False, int(overlap)))
                        temped.write(contigDict[a].forseq + contigDict[b].revseq[int(overlap):] + '\n')
                        edgelist.append((a, b, len(contigDict[a].forseq + contigDict[b].revseq[int(overlap):])))
                    else:
                        contigDict[a].to.append((b, False, overlap))
                        contigDict[b].to.append((a, False, overlap[::-1].translate(transtab)))
                        temped.write(contigDict[a].forseq + overlap + contigDict[b].revseq + '\n')
                        edgelist.append((a, b, len(contigDict[a].forseq + overlap + contigDict[b].revseq)))
            else:
                if dirb == 'True':
                    if overlap.isdigit():
                        contigDict[a].fr.append((b, True, int(overlap)))
                        contigDict[b].fr.append((a, True, int(overlap)))
                        temped.write(contigDict[a].revseq + contigDict[b].forseq[int(overlap):] + '\n')
                        edgelist.append((a, b, len(contigDict[a].revseq + contigDict[b].forseq[int(overlap):])))
                    else:
                        contigDict[a].fr.append((b, True, overlap))
                        contigDict[b].fr.append((a, True, overlap[::-1].translate(transtab)))
                        temped.write(contigDict[a].revseq + overlap + contigDict[b].forseq + '\n')
                        edgelist.append((a, b, len(contigDict[a].revseq + overlap + contigDict[b].forseq)))
                else:
                    if overlap.isdigit():
                        contigDict[a].fr.append((b, False, int(overlap)))
                        contigDict[b].to.append((a, True, int(overlap)))
                        temped.write(contigDict[a].revseq + contigDict[b].revseq[int(overlap):] + '\n')
                        edgelist.append((a, b, len(contigDict[a].revseq + contigDict[b].revseq[int(overlap):])))
                    else:
                        contigDict[a].fr.append((b, False, overlap))
                        contigDict[b].to.append((a, True, overlap[::-1].translate(transtab)))
                        temped.write(contigDict[a].revseq + overlap + contigDict[b].revseq + '\n')
                        edgelist.append((a, b, len(contigDict[a].revseq + overlap + contigDict[b].revseq)))
            count += 1

templg.close()
temped.close()



first = True
reflist = {}
contiglist = {}
for i in sys.argv[2:]:
    tempref = open(i)
    for line in tempref:
        if line.startswith('>'):
            if first:
                first = False
            else:
                reflist[name] = seq
                contiglist[name] = []
            name = line[1:].split()[0]
            seq = ''
        else:
            seq += line.rstrip()
    tempref.close()
reflist[name] = seq
contiglist[name] = []

refout = open('tempref.fa', 'w')
reflen = {}
for i in reflist:
    refout.write('>' + i + '\n' + reflist[i] + reflist[i] + '\n')
    reflen[i] = len(reflist[i])
refout.close()


subprocess.Popen('makeblastdb -dbtype nucl -out tempdb -in tempref.fa', shell=True).wait()
subprocess.Popen('blastn -task blastn -db tempdb -outfmt 6 -query templg.fa -out query_tempdb1.out', shell=True).wait()
subprocess.Popen('blastn -task blastn -db tempdb -outfmt 6 -query tempedge.fa -out query_tempdb2.out', shell=True).wait()

bout = open('query_tempdb1.out')
inref = set()

for line in bout:
    query, subject, ident, length, mismatch, indel, qStart, qEnd, rStart, rEnd, eVal, bitScore = line.split()
    ident = float(ident)
    length = int(length)
    qStart = int(qStart)
    qEnd = int(qEnd)
    rStart = int(rStart)
    rEnd = int(rEnd)
    if ident >= theident and qStart == 1 and qEnd == contigDict[query].length:
        inref.add(query)
        contiglist[subject].append((query, min([rStart, rEnd]), max([rStart, rEnd]), rStart < rEnd, ident, length))
bout.close()
bout = open('query_tempdb2.out')
testset = set()
allconts = set(contigDict)
edges = set()
edgedict ={}

for i in contiglist:
    hitlist = contiglist[i]
    hitlist.sort(key=lambda x: x[-2], reverse=True)
    hitlist.sort(key=lambda x: x[-1], reverse=True)
    newhitlist = []
    for j in hitlist:
        getit = True
        for k in newhitlist:
            if j[1] >= k[1] and j[2] <= k[2]:
                getit = False
                break
        if getit:
            newhitlist.append(j)
    newhitlist.sort(key=lambda x: x[1])
    hitlist = newhitlist
    lasthit = None
    for j in hitlist:
        if lasthit != None:
            if not (j[0], not j[3], lasthit[0], not lasthit[3]) in edges and j[1] < lasthit[2] + 301:
                edges.add((lasthit[0], lasthit[3], j[0], j[3]))
                edgedict[(lasthit[0], lasthit[3], j[0], j[3])] = [lasthit, j]
        lasthit = j

SP = 0
for i in edges:
    #print i
    gotit = False
    if i[1]:
        for j in contigDict[i[0]].to:
            if j[0] == i[2] and j[1] == i[3]:
                gotit = True
                break
    else:
        for j in contigDict[i[0]].fr:
            if j[0] == i[2] and j[1] == i[3]:
                gotit = True
                break
    if gotit:
        SP += 1






edgeinref = set()
for line in bout:
    query, subject, ident, length, mismatch, indel, qStart, qEnd, rStart, rEnd, eVal, bitScore = line.split()
    ident = float(ident)
    length = int(length)
    qStart = int(qStart)
    qEnd = int(qEnd)
    rStart = int(rStart)
    rEnd = int(rEnd)
    if ident >= theident and length >= 0.98 * edgelist[int(query)][-1]:
        edgeinref.add(int(query))



TP = 0
FP = 0
MA = 0
for i in range(len(edgelist)):
    if edgelist[i][0] in inref and edgelist[i][1] in inref:
        if i in edgeinref:
            TP += 1
        else:
            FP += 1
    else:
        MA += 1

print TP , FP , SP, len(edges) - SP
print TP * 100.0 / (TP+FP), '% precision'
print SP * 100.0 / len(edges), '% sensitivity'
