# -*- coding: utf-8 -*-
"""
Created on Dec 21 16:22:14 2019

@author: Zhaoyue Zhang
"""

def printHelpInfo():
    print("""
Useful: 
    python BinomialDistribution.py -k 3 -t DNA -f sequencefile
description:
    extract the number of frequence feature 
    -k k-mer used to extract feature
    -t the type of sequence, can select ['DNA', 'RNA', 'protein']
    -f name of the sequence file folder, storing sequences files with different label
    """)
    exit(0)

def detectingPythonVersionAndOutputHelpInfo():
    if sys.version[0] == '2':
        print("""\nVersionError: The current python version: '%s',
        You must use python version 3 or later to run this script!\n"""%((sys.version).split('\n')[0]))
        exit(0)
    else:
        pass

    try:
        if sys.argv[1] == "--help":
            printHelpInfo()
        else:
            pass
    except:
        printHelpInfo()


def obtainExternalParameters():
    print(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4], sys.argv[5],sys.argv[6])
    try:
        if sys.argv[1] == "-k":
            numConjoin = int(sys.argv[2])
            print(numConjoin)
        else:
            assert 0
        if sys.argv[3] == "-t":
            seqType = str(sys.argv[4])
            print(seqType)
        else:
            assert 0
        if sys.argv[5] == "-f":
            data_path = "%s/"%sys.argv[6]
            print(data_path)
        else:
            assert 0
        
    except:
            printHelpInfo()
    return numConjoin,seqType,data_path

def ObtainKnucleotides(k, seqType):
    if seqType == 'DNA':
        bases = ['A', 'T', 'C', 'G']
        kbases_list_final = allkmer(k, bases)
        return kbases_list_final
    elif seqType == 'RNA':
        bases = ['A', 'U', 'C', 'G']
        kbases_list_final = allkmer(k, bases)
        return kbases_list_final
    elif seqType == 'protein':
        bases = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
        kbases_list_final = allkmer(k, bases)
        return kbases_list_final

def allkmer(k, listseq):
    old = []
    for i in range(k):
        old.append(listseq)
    kbases_list = []
    for pairs in AllPairs(old[:2]):
        kbases_list.append(''.join(pairs))
    for i in range(2, len(old)):
      new_r= []
      for pairs in AllPairs([kbases_list, old[i]]):
        new_r.append(''.join(pairs))
        kbases_list=new_r
    return kbases_list

def ObtainBasesSequence(filename):
    file = open(filename, 'r')
    lines = file.readlines()
    sequences = []
    
    for line in lines:
        if line[0] != '>':
            each_line = line.strip()
            sequences.append(each_line)
    
    return sequences
    
def ObtainSamplesAndLabels(direct_path):
    sample_files = os.listdir(direct_path)
    print(sample_files)
    
    samples_sequences = []
    samples_labels = []
    labels = []
    for files in sample_files:
        seqs = ObtainBasesSequence(direct_path + files)
        samples_sequences.append(seqs)
        seqs_label = len(seqs) * [sample_files.index(files)]
        samples_labels.append(seqs_label)
        labels.append(sample_files.index(files))
    
    return samples_sequences, samples_labels, labels


def Seqtolinesvm(seq, kbases_seq_count, kbases_list, label):
    svmline = str(label)
    
    seq_bases_num = sum(kbases_seq_count.values())
    for kbese in kbases_list:
        num = kbases_seq_count[kbese]
        if num == 0:
            svmline += ',0'
        elif num > 0:
            svmline += ',%.6f'%(num / seq_bases_num)
    
    return svmline + '\n'
        

def CalculateBasesOccur(samples_sequences, labels, kbases_list, featurefile):
    kbases_all_count = dict()
    feature_col = ','.join(kbases_list)
    featurefile.write("label," + feature_col + '\n')    
    
    i_l = 0
    for label in labels:
        kbases_type_count = dict()
        kbases_type_count = kbases_type_count.fromkeys(kbases_list, 0)  
        print("label index: ", labels.index(label))
        samples = samples_sequences[labels.index(label)]
        print("length of samples: ", len(samples))
        j = 0
        for seqs in samples:
            seq = seqs.strip()
            kbases_seq_count = dict()
            kbases_seq_count = kbases_seq_count.fromkeys(kbases_list, 0) 
            kk = 0
            for kbase in kbases_list:
                reg=re.compile("(?=%s)"%(kbase))
                kbase_num = len(reg.findall(seq))
                kbases_seq_count[kbase] = kbase_num
                kbases_type_count[kbase] += kbase_num
                kk = kk + 1
            j = j + 1
            
            featurefile.write(Seqtolinesvm(seq, kbases_seq_count, kbases_list, label))
        kbases_all_count[label] = kbases_type_count
    
    featurefile.close()
    return kbases_all_count
        
def DataskbasesOccur(kbases_all_count, kbases_list, labels):
    datas_kbases = dict()
    datas_kbases = datas_kbases.fromkeys(kbases_list, 0)
    
    for kbase in kbases_list:
        datas_kbases[kbase] = sum([kbases_all_count[label][kbase] for label in labels])
        
    return datas_kbases

def TypePrio(kbases_all_count, datas_kbases, labels):
    prio = dict()
    prio = prio.fromkeys(labels, 0)
    
    M = sum(datas_kbases.values())
    print("Caculating...")
    for label in labels:
        mj = sum(kbases_all_count[label].values())
        prio[label] = mj / M
    
    return prio

def Factorial(n):  
    if n < 0:
        print("Negative numbers have no factorial!")
        return 0
    elif n == 0:
        return 1
    else:
        factorial = 1
        for i in range(1, n + 1):
            factorial = factorial * i
            
        return factorial

def _mainFormula(m, Ni, qj):
    temp = Ni-m
    if Ni != 0:
        if m < temp:
            _nume = range(temp+1, Ni+1)
            _deno = range(1, m+1)
            _grou = zip(_nume, _deno)

            resu = 1

            for (i,j) in _grou:
                resu *= Decimal(((i/j) * qj))
            resu = Decimal(resu * Decimal(((1-qj) ** temp)))
        else:
            _nume = range(m+1, Ni+1)
            _deno = range(1, temp+1)
            _grou = zip(_nume, _deno)

            resu = 1
            for (i,j) in _grou:
                resu *= Decimal(((i/j) * (1-qj)))
            resu = float(Decimal(resu * Decimal((qj ** m))))
            
    else:
        return 1

    return resu

        
def CalculatePnij(kbases_all_count, datas_kbases, labels, kbases_list, prio):
    Pnij = dict()
    ff = open("clresult.txt", "a+")
    
    i_l = 1
    for label in labels:
        qj = prio[label]
        Pni = dict()
        Pni = Pni.fromkeys(kbases_list, 0)
        i_k = 1
        for kbases in kbases_list:
            Ni = datas_kbases[kbases]
            nij = kbases_all_count[label][kbases]
            pni = 0
            for m in range(nij, Ni + 1):
                pni += float(_mainFormula(m, Ni, qj))
            sr = "i_l: " + str(i_l) + "  i_k: " + str(i_k) + "  Ni: " + str(Ni) + "  nij: " + str(nij) + " pnij: " + str(pni) + "\n"
            ff.write(sr)
            i_k = i_k + 1
            Pni[kbases] = pni
        Pnij[label] = Pni
        i_l += 1
    ff.close()
    return Pnij

def CreateCLDict(Pnij, labels, kbases_list):
    CL = dict()
    
    for label in labels:
        CLij = dict()
        CLij = CLij.fromkeys(kbases_list, 0)
        for kbase in kbases_list:
            pnij = Pnij[label][kbase]
            clij = 1 - pnij
            CLij[kbase] = clij
        CL[label] = CLij
    
    return CL

def SortedMaxClvalue(CLDict, labels, kbases_list):
    typelabels = sorted(CLDict.keys())   # 将label升序排列

    kbases_max_CL = dict()
    kbases_max_CL = kbases_max_CL.fromkeys(kbases_list, 0)
    for kbase in kbases_list:
        kbase_max_cl = max([CLDict[label][kbase] for label in typelabels])
        kbases_max_CL[kbase] = kbase_max_cl
        
    return kbases_max_CL

def SaveSortedCL(kbases_max_CL, CLDict, sortedCLfile):
    sorted_cl = sorted(kbases_max_CL.items(), key=lambda asd:asd[1], reverse=True)

    
    filestr = "CLvalue\nfearture\t"
    typelabels = sorted(CLDict.keys())
    for label in typelabels:
        filestr += "label_" + str(label) + '\t'
    filestr += "max_cl\trank\n"
    sortedCLfile.write(filestr)
    
    sorted_kbase = []
    for index_cl in range(len(sorted_cl)):
        (kbase, kbase_cl) = sorted_cl[index_cl]
        sorted_kbase.append(kbase)
        filestr = kbase + "\t"
        for label in typelabels:
            filestr += str("%.6f" %(CLDict[label][kbase])) + "\t"
        filestr += str("%.6f" % (kbases_max_CL[kbase])) + "\t" + str(index_cl) + "\n"
        sortedCLfile.write(filestr)

    sortedCLfile.close()
    
    return sorted_kbase

def getrank(ranked_feature,colnames):
    ranklist = []
    for each in ranked_feature:
        ranklist.append(colnames.index(each))
    return ranklist

def SortedFeaturebyCL(sorted_kbase, featurefilename, sortedfeaturefile):

    f = open(featurefilename, mode = 'r')
    g = open(sortedfeaturefile, mode = 'w')

    countline = 1
    for eachline in f:
        temp_line = eachline.strip().split(',')
        if countline == 1:
            colnames = temp_line[1:]
            ranklist = getrank(sorted_kbase,colnames)
            tempstr = 'label'
            for each in sorted_kbase:
                tempstr +=',%s'%each
            tempstr += '\n'
            countline += 1
            g.write(tempstr)
        else:
            tempstr = '%s'%temp_line[0]
            for each in ranklist:
                tempstr += ',%s'%temp_line[each+1]
            tempstr += '\n'
            g.write(tempstr)

    f.close()
    g.close()
    return


def BinomialFeature(k, data_path, seqType, featurefilename, sortedfeaturefilename, sortedCLfilename):
    featurefile = open(featurefilename, mode = 'w')
    sortedCLfile = open(sortedCLfilename, mode = 'w')
    kbases_list = ObtainKnucleotides(k, seqType)
    samples_sequences, samples_labels, labels = ObtainSamplesAndLabels(data_path)
    kbases_all_count = CalculateBasesOccur(samples_sequences, labels, kbases_list, featurefile)
    datas_kbases = DataskbasesOccur(kbases_all_count, kbases_list, labels)
    prio = TypePrio(kbases_all_count, datas_kbases, labels)
    Pnij = CalculatePnij(kbases_all_count, datas_kbases, labels, kbases_list, prio)
    CLDict = CreateCLDict(Pnij, labels, kbases_list)
    kbases_max_CL = SortedMaxClvalue(CLDict, labels, kbases_list)
    sorted_kbase = SaveSortedCL(kbases_max_CL, CLDict, sortedCLfile)
    SortedFeaturebyCL(sorted_kbase, featurefilename, sortedfeaturefilename)
    print("Finished!")

import sys
import itertools
import re
import pandas as pd
import numpy as np
import os
from decimal import Decimal
import math
from allpairspy import AllPairs

if __name__ == '__main__':
    [numConjoin,seqType,data_path] = obtainExternalParameters()
    BinomialFeature(numConjoin, data_path, seqType, "feature.csv", "BDSortedfeature.csv", "BinomialDistributionCL.txt")

