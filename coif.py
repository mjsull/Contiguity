import networkx
import argparse

class contig:
    def __init__(self):
        pass


def readCSAG(filename):
    pass

# build networkx graph
def buildnx():
    pass

def get_coverage():
    pass

def predict_cov():
    pass

# get length of path
def get_length(path):
    pass

# use all references
def predict_blast1(filename):
    pass

# find best reference
def predict_blast2(filename):
    pass

# read list of predictions
def predict_list(fiename):
    pred = open(filename)
    candidates = {}
    for line in pred:
        candidates[line.split()[0]] = line.split()[1]
    return candidates

# check paths with paired end mapping
def checkPaths(paths):
    pass

# single blast file
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

# scaffold
def scaffold(candContigs):
    pass

def best_guess(candContigs):
    pass

def main():
    readCSAG(filename)
    buildnx()
    if not predictCov is None:
        predFile = predict_cov()
    elif not predictBlast1 is None:
        predFile = predict_blast1(predictBlast1)
    elif not predictBlast2 is None:
        predFile = predict_blast2(predictBlast2)
    elif not predictBlast3 is None:
        predFile = predict_blast3(predictBlast3)
    candContigs = predict_list(predFile)
    if not improvePredict is None:
        predFile = improve_predict(candContigs)
        candContigs = predict_list(predFile)
