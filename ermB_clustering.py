# -*- coding: utf-8 -*-


fi = [line.strip('\n').strip() for line in open('all.align_1-211.list', 'r') ]

# alignment encoding process
# A (1,0,0,0)
# T (0,1,0,0)
# C (0,0,1,0)
# G (0,0,0,1)
# -/ /other (0,0,0,0)
query = []
encodingAln = []
insert = []
seqNum = 0
insertBegin = False
insertDir = {}
rowIndex = -1
for line in fi:
    rowIndex += 1
    if line.startswith('Query_'):
        # define start index
        query_startSeq = line.split()[2][0:10]
        query_len = len(line.split()[2])
        query_startIndex = line.index(query_startSeq)
        # encode query sequence
        for bp in line[query_startIndex:(query_len + query_startIndex)]:
            if bp.upper() == 'A' :
                query.append([1,0,0,0])
            elif bp.upper() == 'T' :
                query.append([0,1,0,0])
            elif bp.upper() == 'C' :
                query.append([0,0,1,0])
            elif bp.upper() == 'G' :
                query.append([0,0,0,1])
            else:
                query.append([0,0,0,0])
        encodingAln.append(query)
    elif line.startswith('Subject'):
        seqNum += 1
        subject = query.copy()
        query_index = -1
        # encode subject sequences but not insert
        for n in range(query_startIndex,(query_len + query_startIndex)):
            bp = line[n]
            query_index += 1
            if bp == '.':
                pass
            elif bp.upper() == 'A' :
                subject[query_index] = [1,0,0,0]
            elif bp.upper() == 'T' :
                subject[query_index] = [0,1,0,0]
            elif bp.upper() == 'C' :
                subject[query_index] = [0,0,1,0]
            elif bp.upper() == 'G' :
                subject[query_index] = [0,0,0,1]
            else:
                subject[query_index] = [0,0,0,0]
        encodingAln.append(subject)
    elif line != "" and line.split()[0].startswith('\\'):
        # process insertion
        insertBegin = True
        insertIndexList = []
        for n in range(query_startIndex, len(fi[rowIndex+2]) ):
            if fi[rowIndex+2][n] == '|':
                insertIndexList.append(n-query_startIndex)
        add = 2
        while insertBegin == True:
            row = rowIndex + add
            insertLine = fi[row].split()
            if len(insertLine) == len(insertIndexList):
                for e in range(len(insertIndexList)):
                    if insertIndexList[e] in insertDir:
                        insertDir[insertIndexList[e]].append((seqNum, insertLine[e]))
                        if len(insertLine[e]) > insertDir[insertIndexList[e]][0]:
                            insertDir[insertIndexList[e]][0] = len(insertLine[e])
                    else:
                        insertDir[insertIndexList[e]] = [len(insertLine[e]), (seqNum, insertLine[e])]
            else:
                insertItem = ''
                for n in range(query_startIndex,len(fi[row])):
                    if  fi[row][n] == "|" or fi[row][n] == " " :
                        if insertItem != "":
                            if (n-1-query_startIndex) in insertDir:
                                insertDir[n-1-query_startIndex].append((seqNum, insertItem))
                                if len(insertItem) > insertDir[n-1-query_startIndex][0]:
                                    insertDir[n-1-query_startIndex][0] = len(insertItem)
                            else:
                                insertDir[n-1-query_startIndex] = [len(insertItem), (seqNum, insertItem)]

                            insertItem = ''
                        else:
                            pass
                    else:
                        insertItem += fi[row][n]
                if insertItem != '':
                    if (len(fi[row])-1-query_startIndex) in insertDir:
                        insertDir[len(fi[row])-1-query_startIndex].append((seqNum, insertItem))
                        if len(insertItem) > insertDir[len(fi[row])-1-query_startIndex][0]:
                            insertDir[len(fi[row])-1-query_startIndex][0] = len(insertItem)
                    else:
                        insertDir[len(fi[row])-1-query_startIndex] = [len(insertItem), (seqNum, insertItem)]
                    insertItem = ''
            add += 2
            if fi[row+1].startswith('Subject') or fi[row+1] == "":
                insertBegin = False
temp = list(insertDir.keys())
temp.sort()
sort_insertDir = {i:insertDir[i] for i in temp}
oldinsertLen = 0

for e,item in sort_insertDir.items():
    insertLen = item[0]
    addList = item[1:].copy()
    num = -1
    addIndex = e + oldinsertLen
    oldinsertLen += insertLen
    for line in encodingAln:
        num += 1
        if len(addList) == 0:
            for q in range(insertLen):
                line.insert(addIndex, [0,0,0,0])
        elif num == addList[0][0]:
            if len(addList[0][1]) == insertLen:
                for bp in addList[0][1][::-1]:
                    if bp.upper() == 'A' :
                        line.insert(addIndex, [1,0,0,0])
                    elif bp.upper() == 'T' :
                        line.insert(addIndex, [0,1,0,0])
                    elif bp.upper() == 'C' :
                        line.insert(addIndex, [0,0,1,0])
                    elif bp.upper() == 'G' :
                        line.insert(addIndex, [0,0,0,1])
                    else:
                        line.insert(addIndex, [0,0,0,0])
            else:
                for q in range(insertLen - len(addList[0][1])):
                    line.insert(addIndex, [0,0,0,0])
                for bp in addList[0][1][::-1]:
                    if bp.upper() == 'A' :
                        line.insert(addIndex, [1,0,0,0])
                    elif bp.upper() == 'T' :
                        line.insert(addIndex, [0,1,0,0])
                    elif bp.upper() == 'C' :
                        line.insert(addIndex, [0,0,1,0])
                    elif bp.upper() == 'G' :
                        line.insert(addIndex, [0,0,0,1])
                    else:
                        line.insert(addIndex, [0,0,0,0])
            addList = addList[1:]
        else:
            for q in range(insertLen):
                line.insert(addIndex, [0,0,0,0])

encodingList = []
for line in encodingAln:
    row = []
    for item in line:
        row.extend(item)
    encodingList.append(row)

import pandas as pd
file = pd.DataFrame(encodingList)

from sklearn.decomposition import PCA
pca = PCA(n_components = 0.95)  # Let M = 2
pca.fit(file)
reduce_data = pca.transform(file)

from sklearn.cluster import KMeans

df = pd.DataFrame()

for n in range(20,101) :
  kmeans = KMeans(n)
  clusters = kmeans.fit_predict(reduce_data)
  df['c'+str(n)] = clusters

df.to_csv('kmeans_211.csv', index=False)