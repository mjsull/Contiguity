import sys, subprocess

theident = 85
bignodesize = 49

class contig:
    def __init__(self, name, forseq, revseq):
        self.name = name
        self.forseq = forseq
        self.revseq = revseq
        self.to = []
        self.fr = []
        self.length = len(forseq)

lg = open(sys.argv[1])
getline = 0
contigDict = {}
edgelist = []
templg = open('templg.fa', 'w')
bignodes = set()
count = 0
for line in lg:
    if line.startswith('NODE'):
        name = line.split()[1]
        getline = 1
    elif getline == 1:
        forseq = line.rstrip()
        getline = 2
    elif getline == 2:
        revseq = line.rstrip()
        getline = 0
        aninstance = contig(name, forseq, revseq)
        if len(forseq) >= bignodesize:
            bignodes.add(name)
            templg.write('>' + name + '\n' + forseq + '\n')
        contigDict[name] = aninstance
    elif line.startswith('ARC'):
        to, fr = line.split()[1:3]
        if to[0] == '-':
            if fr not in contigDict[to[1:]].fr:
                contigDict[to[1:]].fr.append(fr)
        else:
            if fr not in contigDict[to].to:
                contigDict[to].to.append(fr)
        if to[0] == '-':
            to = to[1:]
        else:
            to = '-' + to
        if fr[0] == '-':
            if to not in contigDict[fr[1:]].to:
                contigDict[fr[1:]].to.append(to)
        else:
            if to not in contigDict[fr].fr:
                contigDict[fr].fr.append(to)
templg.close()
temped = open('tempedge.fa', 'w')
todo = []
outpaths = []
for i in bignodes:
    for j in contigDict[i].to:
        todo.append([i, j])
    for j in contigDict[i].fr:
        todo.append(['-' + i, j])

while len(todo) != 0:
    currpath = todo.pop()
    if currpath[-1].replace('-', '') in bignodes:
        outpaths.append(currpath)
    elif len(currpath) < 300:
        if currpath[-1][0] == '-':
            for i in contigDict[currpath[-1][1:]].fr:
                todo.append(currpath + [i])
        else:
            for i in contigDict[currpath[-1]].to:
                todo.append(currpath + [i])

temppaths = []
for i in outpaths:
    if len(i) > 2:
        lall = 0
        for j in i[1:-1]:
            if j[0] == '-':
                lall += contigDict[j[1:]].length
            else:
                lall += contigDict[j].length
        if lall <= 300:
            temppaths.append(i)
        else:
            print 'ding'
    else:
        temppaths.append(i)
       # print i

count = 0
for i in contigDict:
    contigDict[i].to = []
    contigDict[i].fr = []
for i in temppaths:
    seq = ''
    for j in i:
        if j[0] == '-':
            seq += contigDict[j[1:]].revseq
        else:
            seq += contigDict[j].forseq
    edgelist.append((i[0].replace('-', ''), i[-1].replace('-', ''), len(seq)))
    if i[0][0] == '-':
        a = i[0][1:]
        if i[-1][0] == '-':
            b = i[-1][1:]
            contigDict[a].fr.append((b, False))
            contigDict[b].to.append((a, True))
        else:
            b = i[-1]
            contigDict[a].fr.append((b, True))
            contigDict[b].fr.append((a, True))
    else:
        if i[-1][0] == '-':
            b = i[-1][1:]
            contigDict[a].to.append((b, False))
            contigDict[b].to.append((a, False))
        else:
            b = i[-1]
            contigDict[a].to.append((b, True))
            contigDict[b].fr.append((a, False))
    temped.write('>' + str(count) + '\n' + seq + '\n')
    count += 1
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
    else:
        print i
        print edgedict[i]
        print contigDict[i[0]].fr
        print contigDict[i[0]].to






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