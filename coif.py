import networkx
import argparse
import os
import subprocess
import string
import sys
import numpy as np

transtab = string.maketrans('atcgATCG', 'tagcTAGC')


# class containing all information for a contig
class contig:
    def __init__(self, name, shortname, sequence, revseq=None, coverage=None):
        self.name = name
        self.shortname = shortname
        self.forseq = sequence.upper()
        if revseq == None:
            tempseq = self.forseq[::-1]
            self.revseq = tempseq.translate(transtab)
        else:
            self.revseq = revseq
        self.length = len(sequence)
        self.xlength = None
        if self.length >= 1000000000:
            self.strlen = str(round(self.length * 1.0 / 1000000000, 2)) + 'Gb'
        elif self.length >= 1000000:
            self.strlen = str(round(self.length * 1.0 / 1000000, 2)) + 'Mb'
        elif self.length >= 1000:
            self.strlen = str(self.length / 1000) + 'Kb'
        else:
            self.strlen = str(self.length) + 'bp'
        self.visible = False
        self.to = []
        self.fr = []
        if coverage is None:
            self.coverage = 'N/A'
        else:
            self.coverage = round(coverage, 2)
        try:
            self.coverage = float(name.split('_')[5])
        except:
            pass
        gcount = self.forseq.count('G')
        ccount = self.forseq.count('C')
        acount = self.forseq.count('A')
        tcount = self.forseq.count('T')
        self.gccontent = round((gcount + ccount) * 100.0 / self.length, 2)
        try:
            self.gcskew = round((gcount - ccount) * 1.0 / (gcount + ccount), 2)
        except ZeroDivisionError:
            self.gcskew = 0
        try:
            self.atskew = round((acount - tcount) * 1.0 / (acount + tcount), 2)
        except ZeroDivisionError:
            self.atskew = 0


def load_fasta(self):
    fasta = open(args.basic)
    for line in fasta
        if line.startswith('>'):
            if first:
                first = False
            else:
                aninstance = contig(name, name, seq)
                contigDict[name] = aninstance
            name = line.rstrip()[1:]
            seq = ''
        else:
            seq += line.rstrip()
    aninstance = contig(name, name, seq)
    contigDict[name] = aninstance

# load a CAG file
def load_cag(self):
    csag = open(args.cag_file)
    edgelist = []
    global contigDict
    for line in csag:
        if line.split()[0] == 'NODE':
            if len(line.split()) == 4:
                title, entry, name, seq = line.split()
                aninstance = contig(entry, name, seq)
            else:
                title, entry, name, coverage, seq = line.split()
                if coverage == 'N/A':
                    coverage = None
                else:
                    coverage = float(coverage)
                aninstance = contig(entry, name, seq, None, coverage)
            contigDict[entry] = aninstance
        elif line.split()[0] == 'EDGE':
            title, n1, d1, n2, d2, overlap = line.split()
            if overlap == '.':
                overlap = ''
            if d1 == 'True':
                d1 = True
            else:
                d1 = False
            if d2 == 'True':
                d2 = True
            else:
                d2 = False
            if overlap.isdigit():
                overlap = int(overlap)
            edgelist.append((n1, d1, n2, d2, overlap))
    for i in edgelist:
        contiga, dira, contigb, dirb, overlap = i
        if dira and dirb:
            contigDict[contiga].to.append((contigb, True, overlap))
            contigDict[contigb].fr.append((contiga, False, overlap))
        elif dira and not dirb:
            contigDict[contiga].to.append((contigb, False, overlap))
            contigDict[contigb].to.append((contiga, False, overlap))
        elif not dira and dirb:
            contigDict[contiga].fr.append((contigb, True, overlap))
            contigDict[contigb].fr.append((contiga, True, overlap))
        else:
            contigDict[contiga].fr.append((contigb, False, overlap))
            contigDict[contigb].to.append((contiga, True, overlap))

# find best reference
def predict_blast1(args):
    maxgot = 0
    maxgot2 = 0
    bestref2 = None
    bestref = None
    for i in os.listdir(args.blast1):
        reflen = 0
        ref = open(args.blast1 + '/' + i)
        for line in ref:
            if not line.startswith('>'):
                reflen += len(line.rstrip())
        ref.close()
        subprocess.Popen('makeblastdb -dbtype nucl -out ' + args.work_dir + '/tempdb -in ' + args.blast1 + '/' + i, shell=True, stdout=subprocess.PIPE).wait()
        subprocess.Popen('blastn -db ' + args.work_dir + '/tempdb -outfmt 6 -num_threads 8 -query ' + args.work_dir + '/contigs.fa -out ' + args.work_dir + '/contigs_tempdb.out', shell=True).wait()
        blast = open(args.work_dir + '/contigs_tempdb.out')
        gotset = set()
        querydict = {}
        for line in blast:
            query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = line.split()
            qstart, qstop, rstart, rstop, length, mm = map(int, [qstart, qstop, rstart, rstop, length, mm])
            eval = float(eval)
            if eval <= 0.005:
                for j in range(min([rstart, rstop]), max([rstart, rstop]) + 1):
                    gotset.add(j)
                if not query in querydict:
                    querydict[query] = set()
                for j in range(qstart, qstop + 1):
                    querydict[query].add(j)
        blast.close()
        aval, bval = 0, 0
        for j in contigDict:
            if j in querydict:
                aval += len(querydict[j])
            bval += contigDict[j].length
        gotset2 = aval * 1.0 / bval
        if len(gotset) * 1.0/ reflen + gotset2 >= maxgot:
            bestref2 = bestref
            maxgot2 = maxgot
            bestref = i
            maxgot = len(gotset) * 1.0 / reflen + gotset2
        elif len(gotset) * 1.0/reflen + gotset2 >= maxgot2:
            bestref2 = i
            maxgot2 = len(gotset) * 1.0 /reflen + gotset2
    return bestref2

# write CAG to FASTA for perfoming BLAST etc.
def write_fasta_cag(args):
    out = open(args.work_dir + '/contigs.fa', 'w')
    for i in contigDict:
        out.write('>' + i + '\n' + contigDict[i].forseq + '\n')
    out.close()

# use selected reference to predict plasmid contigs
def predict_blast2(args):
    subprocess.Popen('makeblastdb -dbtype nucl -out ' + args.work_dir + '/tempdb -in ' + args.blast2, shell=True, stdout=subprocess.PIPE).wait()
    subprocess.Popen('blastn -db ' + args.work_dir + '/tempdb -outfmt 6 -num_threads 8 -query ' + args.work_dir + '/contigs.fa -out ' + args.work_dir + '/contigs_tempdb.out', shell=True).wait()
    filtered = set()
    blast = open(args.work_dir + '/contigs_tempdb.out')
    for line in blast:
        query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = line.split()
        qstart, qstop, rstart, rstop, length, mm = map(int, [qstart, qstop, rstart, rstop, length, mm])
        ident = float(ident)
        if length >= args.min_length and ident >= args.min_ident and length >= contigDict[query].length * args.min_len_fract:
            filtered.add(query)
    blast.close()
    candidates = set()
    for i in contigDict:
        if not i in filtered:
            candidates.add(i)
    return candidates

# turn the contig dictionary into a Directional graph for networkx
def contigDict2nx():
    dg = networkx.DiGraph()
    for i in contigDict:
        dg.add_node(i)
        dg.add_node('-' + i)
        for j in contigDict[i].to:
            if j[1]:
                dg.add_edge(i, j[0])
            else:
                dg.add_edge(i, '-' + j[0])
        for j in contigDict[i].fr:
            if j[1]:
                dg.add_edge('-' + i, j[0])
            else:

                dg.add_edge('-' + i, '-' + j[0])
    return dg

# Find all simple paths between a list of candidates
def find_paths(dg, candidates, args):
    outpaths = []
    listcan = list(candidates)
    while len(listcan) > 0:
        i = listcan[0]
        for j in listcan:
            paths = list(networkx.all_simple_paths(dg, i, j, args.max_path_node))
            for k in paths:
                pathseq = ''
                lastcontig = i
                lastcontigdir = True
                getit = True
                for l in k[1:-1]:
                    if '-' in l:
                        contig = l[1:]
                        contigdir = False
                    else:
                        contig = l
                        contigdir = True
                    if contig in candidates:
                        getit = False
                        break
                    if lastcontigdir:
                        for m in contigDict[lastcontig].to:
                            if m[0] == contig and m[1] == contigdir:
                                overlap = m[2]
                                break
                    else:
                        for m in contigDict[lastcontig].fr:
                            if m[0] == contig and m[1] == contigdir:
                                overlap = m[2]
                                break
                    if type(overlap) is int:
                        if contigdir:
                            pathseq += contigDict[contig].forseq[overlap:]
                        else:
                            pathseq += contigDict[contig].revseq[overlap:]
                    else:
                        pathseq += overlap
                        if contigdir:
                            pathseq += contigDict[contig].forseq
                        else:
                            pathseq += contigDict[contig].revseq
                    lastcontig = contig
                    lastcontigdir = contigdir
                if k[-1][0] == '-':
                    contig = k[-1][1:]
                    contigdir = False
                else:
                    contig = k[-1]
                    contigdir = True
                if lastcontigdir:
                    for m in contigDict[lastcontig].to:
                        if m[0] == contig and m[1] == contigdir:
                            overlap = m[2]
                            break
                else:
                    for m in contigDict[lastcontig].fr:
                        if m[0] == contig and m[1] == contigdir:
                            overlap = m[2]
                            break
                if getit:
                    if type(overlap) is int:
                        pathseq = pathseq[:-overlap]
                    else:
                        pathseq += overlap
                    if len(pathseq) <= args.max_path_length:
                        outpaths.append((k, pathseq))
            paths = list(networkx.all_simple_paths(dg, i, '-' + j, args.max_path_node))
            for k in paths:
                pathseq = ''
                lastcontig = i
                lastcontigdir = True
                getit = True
                for l in k[1:-1]:
                    if '-' in l:
                        contig = l[1:]
                        contigdir = False
                    else:
                        contig = l
                        contigdir = True
                    if contig in candidates:
                        getit = False
                        break
                    if lastcontigdir:
                        for m in contigDict[lastcontig].to:
                            if m[0] == contig and m[1] == contigdir:
                                overlap = m[2]
                                break
                    else:
                        for m in contigDict[lastcontig].fr:
                            if m[0] == contig and m[1] == contigdir:
                                overlap = m[2]
                                break
                    if type(overlap) is int:
                        if contigdir:
                            pathseq += contigDict[contig].forseq[overlap:]
                        else:
                            pathseq += contigDict[contig].revseq[overlap:]
                    else:
                        pathseq += overlap
                        if contigdir:
                            pathseq += contigDict[contig].forseq
                        else:
                            pathseq += contigDict[contig].revseq
                    lastcontig = contig
                    lastcontigdir = contigdir
                if k[-1][0] == '-':
                    contig = k[-1][1:]
                    contigdir = False
                else:
                    contig = k[-1]
                    contigdir = True
                if lastcontigdir:
                    for m in contigDict[lastcontig].to:
                        if m[0] == contig and m[1] == contigdir:
                            overlap = m[2]
                            break
                else:
                    for m in contigDict[lastcontig].fr:
                        if m[0] == contig and m[1] == contigdir:
                            overlap = m[2]
                            break
                if getit:
                    if type(overlap) is int:
                        pathseq = pathseq[:-overlap]
                    else:
                        pathseq += overlap
                    if len(pathseq) <= args.max_path_length:
                        outpaths.append((k, pathseq))
            paths = list(networkx.all_simple_paths(dg, '-' + i, j, args.max_path_node))
            for k in paths:
                pathseq = ''
                lastcontig = i
                lastcontigdir = False
                getit = True
                for l in k[1:-1]:
                    if '-' in l:
                        contig = l[1:]
                        contigdir = False
                    else:
                        contig = l
                        contigdir = True
                    if contig in candidates:
                        getit = False
                        break
                    if lastcontigdir:
                        for m in contigDict[lastcontig].to:
                            if m[0] == contig and m[1] == contigdir:
                                overlap = m[2]
                                break
                    else:
                        for m in contigDict[lastcontig].fr:
                            if m[0] == contig and m[1] == contigdir:
                                overlap = m[2]
                                break
                    if type(overlap) is int:
                        if contigdir:
                            pathseq += contigDict[contig].forseq[overlap:]
                        else:
                            pathseq += contigDict[contig].revseq[overlap:]
                    else:
                        pathseq += overlap
                        if contigdir:
                            pathseq += contigDict[contig].forseq
                        else:
                            pathseq += contigDict[contig].revseq
                    lastcontig = contig
                    lastcontigdir = contigdir
                if k[-1][0] == '-':
                    contig = k[-1][1:]
                    contigdir = False
                else:
                    contig = k[-1]
                    contigdir = True
                if lastcontigdir:
                    for m in contigDict[lastcontig].to:
                        if m[0] == contig and m[1] == contigdir:
                            overlap = m[2]
                            break
                else:
                    for m in contigDict[lastcontig].fr:
                        if m[0] == contig and m[1] == contigdir:
                            overlap = m[2]
                            break
                if getit:
                    if type(overlap) is int:
                        pathseq = pathseq[:-overlap]
                    else:
                        pathseq += overlap
                    if len(pathseq) <= args.max_path_length:
                        outpaths.append((k, pathseq))
            paths = list(networkx.all_simple_paths(dg, '-' + i, '-' + j, args.max_path_node))
            for k in paths:
                pathseq = ''
                lastcontig = i
                lastcontigdir = False
                getit = True
                for l in k[1:-1]:
                    if '-' in l:
                        contig = l[1:]
                        contigdir = False
                    else:
                        contig = l
                        contigdir = True
                    if contig in candidates:
                        getit = False
                        break
                    if lastcontigdir:
                        for m in contigDict[lastcontig].to:
                            if m[0] == contig and m[1] == contigdir:
                                overlap = m[2]
                                break
                    else:
                        for m in contigDict[lastcontig].fr:
                            if m[0] == contig and m[1] == contigdir:
                                overlap = m[2]
                                break
                    if type(overlap) is int:
                        if contigdir:
                            pathseq += contigDict[contig].forseq[overlap:]
                        else:
                            pathseq += contigDict[contig].revseq[overlap:]
                    else:
                        pathseq += overlap
                        if contigdir:
                            pathseq += contigDict[contig].forseq
                        else:
                            pathseq += contigDict[contig].revseq
                    lastcontig = contig
                    lastcontigdir = contigdir
                if k[-1][0] == '-':
                    contig = k[-1][1:]
                    contigdir = False
                else:
                    contig = k[-1]
                    contigdir = True
                if lastcontigdir:
                    for m in contigDict[lastcontig].to:
                        if m[0] == contig and m[1] == contigdir:
                            overlap = m[2]
                            break
                else:
                    for m in contigDict[lastcontig].fr:
                        if m[0] == contig and m[1] == contigdir:
                            overlap = m[2]
                            break
                if getit:
                    if type(overlap) is int:
                        pathseq = pathseq[:-overlap]
                    else:
                        pathseq += overlap
                    if len(pathseq) <= args.max_path_length:
                        outpaths.append((k, pathseq))
        listcan.pop(0)
    return outpaths



# check paths with paired end mapping
def check_paths(paths):
    pass


# given a list of paths find the shortest
def getShortest(paths):
    outpaths = {}
    for i in paths:
        sn = i[0][0]
        en = i[0][-1]
        if sn in outpaths:
            if en in outpaths[sn]:
                if len(i[1]) < len(outpaths[sn][en][1]):
                    outpaths[sn][en] = i
            else:
                outpaths[sn][en] = i
        else:
            outpaths[sn] = {en:i}
    return outpaths



def predict_blast3(filename):
    pass

# improve predictions
def improve_predict(candContigs):
    pass

# predict small contigs
def predict_small(candContigs):
    pass

# remove similar paths
def remove_dup(paths):
    pass


# find all simple circuits using a graph of candidates
def find_circuits(canddict):
    nx = networkx.DiGraph()
    for i in canddict:
        for j in canddict[i]:
            nx.add_edge(i, j)
            if i[0] == '-':
                newi = i[1:]
            else:
                newi = '-' + i
            if j[0] == '-':
                newj = j[1:]
            else:
                newj = '-' + j
            nx.add_edge(newj, newi)
    return networkx.simple_cycles(nx)



def predict_contigs(canddict, plascan):
    sys.stdout.write('Trimming contigs not in circuit..\n')
    G = networkx.DiGraph()
    for i in canddict:
        for j in canddict[i]:
            G.add_edge(i, j)
            if i[0] == '-':
                newi = i[1:]
            else:
                newi = '-' + i
            if j[0] == '-':
                newj = j[1:]
            else:
                newj = '-' + j
            G.add_edge(newj, newi)
    sys.stdout.write(str(len(G.nodes())) + ' initial nodes in ' + str(networkx.number_connected_components(G.to_undirected())) + ' connected components.\n')
    notgotemall = True
    while notgotemall:
        notgotemall = False
        for i in plascan:
            if i in G:
                if len(G[i]) == 0 or len(G['-' + i]) == 0:
                    notgotemall = True
                    G.remove_node(i)
                    G.remove_node('-' + i)
    sys.stdout.write(str(len(G.nodes())) + ' remaining nodes.\n')
    predset = set()
    for i in G.edges():
        if i[0] in canddict and i[1] in canddict[i[0]]:
            for j in canddict[i[0]][i[1]][0]:
                if j[0] == '-':
                    predset.add(j[1:])
                else:
                    predset.add(j)
    return predset

def theil_sen(x,y):
    n = len(x)
    ord = np.argsort(x)
    xs = x[ord]
    ys = y[ord]
    vec1 = np.zeros( (n,n) )
    for ii in range(n):
        for jj in range(n):
            vec1[ii,jj] = ys[ii]-ys[jj]
    vec2 = np.zeros( (n,n) )
    for ii in range(n):
        for jj in range(n):
            vec2[ii,jj] = xs[ii]-xs[jj]
    v1 = vec1[vec2>0]
    v2 = vec2[vec2>0]
    slope = np.median( v1/v2 )
    coef = np.zeros( (2,1) )
    b_0 = np.median(y)-slope*np.median(x)
    b_1 = slope
    res = y-b_1*x-b_0 # residuals
    return (b_0,b_1)


def predict_cov(args):
    plascan = set()
    count = 0
    for i in contigDict:
        if contigDict[i].length >= args.min_can_length:
            count += 1
    x = np.zeros(count)
    y = np.zeros(count)
    contignamelist = []
    index = 0
    for i in contigDict:
        if contigDict[i].length >= args.min_can_length:
            x[index] = contigDict[i].gccontent
            y[index] = contigDict[i].coverage
            contignamelist.append(i)
            index += 1
    thec, ther = theil_sen(x, y)
    lessthan = []
    for i in range(len(x)):
        if x[i] * ther + thec - y[i] >= 0:
            lessthan.append(x[i] * ther + thec - y[i])
    lessthansd = (sum(lessthan) * 1.0 / len(lessthan)) ** 0.5
    templessthan = []
    for i in lessthan:
        if i <= 4 * lessthansd:
            templessthan.append(i)
    lessthan = templessthan
    lessthan.sort()
    thecutoff = lessthan[int(0.95 * len(lessthan))]
    out = open(args.work_dir + '/coverage.csv', 'w')
    for i in contigDict:
        out.write(i + '\t' + str(contigDict[i].coverage) + '\t' + str(contigDict[i].gccontent) + '\t' + str(contigDict[i].length)
                  + '\t' + str(contigDict[i].gccontent * ther + thec + thecutoff) + '\t' + str(contigDict[i].gccontent * ther + thec) + '\n')
    out.write('\n\n\n')
    for i in contigDict:
        if contigDict[i].length > args.min_can_length:
            out.write(i + '\t' + str(contigDict[i].coverage) + '\t' + str(contigDict[i].gccontent) + '\t' + str(contigDict[i].length)
                  + '\t' + str(contigDict[i].gccontent * ther + thec + thecutoff) + '\t' + str(contigDict[i].gccontent * ther + thec) + '\n')
    for i in contigDict:
        contigDict[i].predScore = min([2, (contigDict[i].coverage - (contigDict[i].gccontent * ther + thec)) / thecutoff])
        if contigDict[i].length >= args.min_can_length and contigDict[i].coverage >= contigDict[i].gccontent * ther + thec + thecutoff:
            if args.filter_high_cov and contigDict[i].coverage <= contigDict[i].gccontent * ther * 2 + thec * 2 - thecutoff:
                plascan.add(i)
            elif not args.filter_high_cov:
                plascan.add(i)
    out.close()
    return plascan




# scaffold
def scaffold(candContigs):
    pass

def best_guess(candContigs):
    pass


def get_pred_qual(args, candidates, initit, onlyplas=False):
    first = True
    plasconts = set()
    chromconts = set()
    totalplaslengths = 0
    sharedbp, TPbp, FPbp, TNbp, FNbp, unmappedbp = 0, 0, 0, 0, 0, 0
    for i in args.debug:
        tempref = open(args.work_dir + '/tempseq.fa', 'w')
        inref = open(i)
        for line in inref:
            if line.startswith('>'):
                tempref.write(line)
                seq = ''
            else:
                seq += line.rstrip()
        tempref.write(seq + seq + '\n')
        tempref.close()
        thelen = len(seq)
        gotthem = set()
        gotthem2 = set()
        subprocess.Popen('makeblastdb -dbtype nucl -out ' + args.work_dir + '/tempdb -in ' + args.work_dir + '/tempseq.fa', shell=True, stdout=subprocess.PIPE).wait()
        subprocess.Popen('blastn -db ' + args.work_dir + '/tempdb -outfmt 6 -num_threads 8 -query ' + args.work_dir + '/contigs.fa -out ' + args.work_dir + '/contigs_tempdb.out', shell=True).wait()
        blast = open(args.work_dir + '/contigs_tempdb.out')
        for line in blast:
            query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = line.split()
            qstart, qstop, rstart, rstop, length, mm = map(int, [qstart, qstop, rstart, rstop, length, mm])
            ident = float(ident)
            if length >= 0.99 * contigDict[query].length and ident >= 90.0:
                if first and not onlyplas:
                    chromconts.add(query)
                else:
                    plasconts.add(query)
                    if min([rstart, rstop]) <= thelen:
                        if query in candidates:
                            for q in range(min([rstart, rstop]), max([rstart, rstop]) + 1):
                                gotthem.add(q)
                        for q in range(min([rstart, rstop]), max([rstart, rstop]) + 1):
                            gotthem2.add(q)
        if not first or onlyplas:
            TPbp += len(gotthem)
            FNbp += len(gotthem2)
            totalplaslengths += thelen
        first = False
    shared, TP, FP, TN, FN, unmapped = 0, 0, 0, 0, 0, 0
    for i in contigDict:
        if i in chromconts and i in plasconts:
            shared += 1
            sharedbp += contigDict[i].length
            if i in candidates:
                TP += 1
            else:
                FN += 1
        elif i in chromconts:
            if i in candidates:
                FP += 1
                FPbp += contigDict[i].length
            else:
                TN += 1
                TNbp += contigDict[i].length
        elif i in plasconts:
            if i in candidates:
                TP += 1
            else:
                FN += 1
        else:
            unmapped += 1
            unmappedbp += contigDict[i].length
            if onlyplas:
                if i in candidates:
                    FP += 1
                    FPbp += contigDict[i].length
                else:
                    TN += 1
                    TNbp += contigDict[i].length
    if initit:
        out = open(args.work_dir + '/degbug.txt', 'w')
    else:
        out = open(args.work_dir + '/degbug.txt', 'a')
    #out.write('candidates\n' + '\t'.join(candidates) + '\n')
    if initit:
        out.write('predictive power initial set\ntp\tfp\ttn\tfn\tsensitivity\tprecision\tunmapped\tshared\n')
    else:
        out.write('predictive power final set\ntp\tfp\ttn\tfn\tsensitivity\tprecision\tunmapped\tshared\n')
    try:
        out.write('\t'.join(map(str, [TP, FP, TN, FN, TP * 1.0 / (TP+FN), TP * 1.0 / (TP+FP), unmapped, shared])) + '\n')
    except:
        out.write('\t'.join(map(str, [TP, FP, TN, FN, 0, 0, unmapped, shared])) + '\n')
    if initit:
        out.write('predictive power initial set (bp)\ntp\tfp\ttn\tfn\tsensitivity\tprecision\tunmapped\tshared\tplasass\n')
    else:
        out.write('predictive power final set (bp)\ntp\tfp\ttn\tfn\tsensitivity\tprecision\tunmapped\tshared\tplasass\n')
    try:
        out.write('\t'.join(map(str, [TPbp, FPbp, TNbp, FNbp - TPbp, TPbp * 1.0 / FNbp, TPbp * 1.0 / (TPbp+FPbp), unmappedbp, sharedbp, FNbp * 1.0 / totalplaslengths])) + '\n')
    except:
        out.write('\t'.join(map(str, [TPbp, FPbp, TNbp, FNbp - TPbp, 0, 0, unmappedbp, sharedbp, 0])) + '\n')
    out.close()

# write list of candidates to FASTA file
def write_cand(candidates, outfile):
    out = open(args.work_dir + '/' + outfile, 'w')
    for i in candidates:
        out.write(i + '\n')
    out.close()


def main(args):
    global contigDict
    contigDict = {}
    if args.basic is None:
        load_cag(args)
    else:
        load_fasta(args)
    if os.path.exists(args.work_dir): # create a working directory if it doesn't exist
        if not os.path.isdir(args.work_dir):
            sys.stderr.write('Working directory is a file not a folder.\n')
            sys.exit()
    else:
        os.makedirs(args.work_dir)
    write_fasta_cag(args) # create a FASTA file of the graph file in the working directory
    if args.predict_cov: # if flag set predict initial candidates using coverage
        sys.stdout.write('Finding candidates.\n')
        candidates = predict_cov(args)
        sys.stdout.write(str(len(candidates)) + ' candidates found.\n')
    elif not args.blast1 is None:
        sys.stdout.write('Finding best reference..\n')
        bestref = predict_blast1(args)
        sys.stdout.write('Using reference ' + bestref + '\n')
        args.blast2 = args.blast1 + '/' + bestref
        sys.stdout.write('Finding candidates.\n')
        candidates = predict_blast2(args)
        sys.stdout.write(str(len(candidates)) + ' candidates found.\n')
    elif not args.blast2 is None:
        sys.stdout.write('Finding candidates.\n')
        candidates = predict_blast2(args)
        sys.stdout.write(str(len(candidates)) + ' candidates found.\n')
    if not args.debug is None:
        get_pred_qual(args, candidates, True, args.only_plas)
    write_cand(candidates, 'candidates.txt')
    if args.basic is None:
        sys.stdout.write('Basic mode finished.')
        return
    sys.stdout.write('Creating graph\n')
    dg = contigDict2nx()
    sys.stdout.write('Graph created.\nFinding paths between candidates\n')
    paths = find_paths(dg, candidates, args)
    sys.stdout.write(str(len(paths)) + ' paths found.\n')
    if args.check_paths: # check paths with paired-end reads TODO
        pass
    sys.stdout.write('Finding shortest paths.\n')
    pathDict = getShortest(paths)
    count = 0
    for i in pathDict:
        count += len(pathDict[i])
    sys.stdout.write(str(count) + ' paths remaining.\nPredicting final set of contigs..\n')
    newcandidates = predict_contigs(pathDict, candidates)
    sys.stdout.write(str(len(newcandidates)) + ' predicted plasmid contigs.\n')
    if not args.debug is None:
        get_pred_qual(args, newcandidates, False, args.only_plas)
    sys.stdout.write('Writing predicted contigs to file, thanks for using COIF.\n')
    write_cand(newcandidates, 'final_candidates.txt')






parser = argparse.ArgumentParser(prog='coif.py', formatter_class=argparse.RawDescriptionHelpFormatter, description='''
coif.py: A script for identifying plasmid contigs.

USAGE: coif.py -c assembly.cag -d working_dir -b1 folder_of_reference_genomes


''', epilog="Thanks for using Contiguity")
parser.add_argument('-c', '--cag_file', action='store', help='CAG file of assembled contigs or scaffolds and graph')
parser.add_argument('-b', '--basic', action='store', default=None, help='Only do initial prediction, do not improve predictions with graph information. Uses a FASTA file as input instead of a CAG.')
parser.add_argument('-d', '--work_dir', action='store', help='Working directory')
parser.add_argument('-b1', '--blast1', action='store', default=None, help='Find best reference from folder of references and use it to predict initial set of plasmid contigs.')
parser.add_argument('-b2', '--blast2', action='store', default=None, help='Use reference to remove chromsomal contigs.')
parser.add_argument('-pc', '--predict_cov', action='store_true', default=False, help='Find best reference from folder of references and use it to predict initial set of plasmid contigs.')
parser.add_argument('-db', '--debug', action='store', default=None, nargs='+', help='Give references to report performance of COIF [chromosome plas1 plas2 etc.].')
parser.add_argument('-i', '--min_ident', action='store', type=float, default=80.0, help='Min idenity of hits to draw')
parser.add_argument('-l', '--min_length', action='store', type=int, default=0, help='Min length of hits to draw')
parser.add_argument('-f', '--min_len_fract', action='store', type=float, default=0.1, help='Min length of hits to draw')
parser.add_argument('-mp', '--max_path_length', action='store', type=int, default=15000, help='Max length (bp) of paths')
parser.add_argument('-mn', '--max_path_node', action='store', type=int, default=10, help='Max nodes to search paths')
parser.add_argument('-cp', '--check_paths', action='store_true', default=False, help='Check paths with paired end reads')
parser.add_argument('-op', '--only_plas', action='store_true', default=False, help='only use plasmids for debug mode')
parser.add_argument('-mc', '--min_can_length', action='store', type=int, default=500, help='minimum length of contig for initial predictions')
parser.add_argument('-fh', '--filter_high_cov', action='store_true', default=True, help='Check paths with paired end reads')
# parser.add_argument('-rf', '--read_file', action='store', help='read file')
# parser.add_argument('-o', '--output_folder', action='store', help='output folder')
# parser.add_argument('-k', '--kmer_size', action='store', type=int, default=31, help='k-mer size for finding adjacent contigs [31]')
# parser.add_argument('-max_d', '--max_distance', action='store', type=int, default=300, help='maximum distance apart in the de bruijn graph for contigs to count as adjacent [300]')
# parser.add_argument('-kmer_a', '--kmer_average', action='store', type=int, default=-1, help='All k-mers above half this value will be traversed [auto]')
# parser.add_argument('-kmer_c', '--kmer_cutoff', action='store', type=int, default=-1, help='cutoff for k-mer values [auto]')
# parser.add_argument('-ov', '--overlap', action='store', type=int, default=None, help='minimum overlap to create edge [kmer_size-1]')
# parser.add_argument('-rl', '--min_read_length', action='store', type=int, default=75, help='Minimum read length [75]')
# parser.add_argument('-max_mm', '--max_mismatch', action='store', type=int, default=2, help='maximum number of mismatches to count overlap [2]')
# parser.add_argument('-lo', '--long_overlap_ident', action='store', type=int, default=85, help='minimum percent identity to create an edge where there is a long overlap [85]')
# parser.add_argument('-mp', '--minimum_pairs_edge', action='store', type=int, default=2, help='Minimum pairs to create edge [2]')
# parser.add_argument('-is', '--max_insert_size', action='store', type=int, default=600, help='Upper bound on insert size [600]')
# parser.add_argument('-cl', '--command_line', action='store_true', default=False, help='Run contiguity in command line mode')
# parser.add_argument('-no', '--no_overlap_edges', action='store_true', default=False, help='Don\'t get overlap edges')
# parser.add_argument('-nd', '--no_db_edges', action='store_true', default=False, help='Don\'t get De Bruijn edges')
# parser.add_argument('-np', '--no_paired_edges', action='store_true', default=False, help='Don\'t get paired-end edges')
# parser.add_argument('-km', '--khmer', action='store_false', default=True, help='Don\'t use khmer for De Bruijn graph contruction (not recommended)')
# parser.add_argument('-nt', '--num_threads', action='store', type=int, default=1, help='Number of threads to use for hash table building with khmer and for mapping reads with bowtie')
# parser.add_argument('-ht_s', '--ht_size', action='store', default='2e9', help='Hash table size, for more information check http://khmer.readthedocs.org/en/v1.1/choosing-table-sizes.html')
# parser.add_argument('-ht_n', '--ht_number', action='store', type=int, default=4, help='Hash table number, for more information check http://khmer.readthedocs.org/en/v1.1/choosing-table-sizes.html')



args = parser.parse_args()

main(args)