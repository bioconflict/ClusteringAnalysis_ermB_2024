# comparison calculation
# version 1.2
# only could count mismatch, deletion, and insertion
# Yongjun Tan
# May 28, 2021
fi = open('ermB.align.list','r')
fw = open('ermB.align.table','w')
fileEnd = False
readStart = False
subjectNum = 0
nucl = True
pro = False
while fileEnd == False:
    line = fi.readline().strip('\n')
    if line.startswith('Query_1'):
        readStart = True
        seq = list(line.split()[2].upper())
        seqIndex = [[0,0,{}] for n in range(len(seq))]
        insertIndex = [[0,{}] for n in range(len(seq))]
        startNum = line.index(''.join(seq[0:5]))
    elif readStart == True and len(line) == 0:
        fileEnd = True
    elif readStart == True and not line.startswith(' '):
        subjectNum += 1
        match = list(line)
        for aa in range(startNum,len(match)):
            if match[aa].upper() in ['B','J','O','X','Z',' ',
                                     '0','1','2','3','4','5','6',
                                     '7','8','9']:
                pass
            elif match[aa] == '.':
                position = aa - startNum
                seqIndex[position][1] += 1
            elif match[aa] in ['A','G','T','C','-']:
                position = aa - startNum
                seqIndex[position][0] += 1
                seqIndex[position][1] += 1
                if match[aa].upper() in seqIndex[position][2]:
                    seqIndex[position][2][match[aa].upper()] += 1
                else:
                    seqIndex[position][2][match[aa].upper()] = 1
    elif readStart == True and line.startswith(' ') and line.split()[0][0] == '\\':
        insertPosition = []
        match = list(line)
        for aa in range(startNum,len(match)):
            if match[aa].upper() in ['.','B','J','O','X','Z',' ',
                                     '0','1','2','3','4','5','6',
                                     '7','8','9']:
                pass
            elif match[aa] == '\\' :
                insertPosition.append(aa)
                insertIndex[aa-startNum][0] += 1
        loop = True
        while loop == True:
            loop = False
            line_insert = fi.readline()
            line_insert = fi.readline().strip('\n')
            insertStart = startNum
            
            for aa in insertPosition:
                if aa < len(line_insert) and line_insert[aa] not in['|',' ']:
                    if line_insert[insertStart:aa+1].strip() in insertIndex[aa-startNum][1]:
                        insertIndex[aa-startNum][1][line_insert[insertStart:aa+1].strip()] += 1
                    else:
                        insertIndex[aa-startNum][1][line_insert[insertStart:aa+1].strip()] = 1
                    insertStart = aa +1
                elif aa > len(line_insert):
                    pass
                elif line_insert[aa] == '|' :
                    insertStart = aa +1
                    loop = True
            

for n in range(len(seqIndex)):   
    fw.write(seq[n].upper()+'\t'+str(n+1)+'\t'+str(seqIndex[n][0]/seqIndex[n][1])+'\t'+str(seqIndex[n][1])+'\t'+str(seqIndex[n][0]/subjectNum)+'\t'+str(subjectNum))
    if len(seqIndex[n][2]) > 0:
        for k,v in seqIndex[n][2].items():
            fw.write('\t'+k+':'+str(v))
    fw.write('\n')
fw.write('\n')
for n in range(len(insertIndex)):
    if insertIndex[n][0] >0:
        fw.write(str(n)+':'+str(n+1)+'\t'+str(insertIndex[n][0]/subjectNum))
        for k,v in insertIndex[n][1].items():
            fw.write('\t'+k+':'+str(v))
        fw.write('\n')
fi.close()
fw.close()
        
    
