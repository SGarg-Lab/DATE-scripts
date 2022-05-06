import numpy as np
import sys

separation = 1000

def kmerge(reps):
    result = []
    while(reps!=[]):
        incompleteEnhancer = True
        initializedEnhancer = False
        while (incompleteEnhancer):
            while ([] in reps):
                reps.remove([])

            smallest = float('inf')
            minCheck = []
            for i in range(0,len(reps)):
                if int(reps[i][0][1])<smallest:
                    smallest = int(reps[i][0][1])
                    minCheck = []
                    minCheck.append(i)
                elif int(reps[i][0][1])==smallest:
                    minCheck.append(i)

            if len(minCheck)==0:
                incompleteEnhancer = False
                result.append(thisEnhancer)

            elif len(minCheck)==1:
                if (initializedEnhancer == False):
                    thisEnhancer = reps[minCheck[0]][0]
                    del reps[minCheck[0]][0]
                    initializedEnhancer = True
                elif (int(reps[minCheck[0]][0][1])-int(thisEnhancer[2])>separation):
                    incompleteEnhancer = False
                    result.append(thisEnhancer)
                else:
                    if int(thisEnhancer[2]) < int(reps[minCheck[0]][0][2]):
                        thisEnhancer[2] = reps[minCheck[0]][0][2]
                    if thisEnhancer[3] != reps[minCheck[0]][0][3]:
                        thisEnhancer[3] = 'mixed - ERROR'
                    if thisEnhancer[4] != reps[minCheck[0]][0][4]:
                        thisEnhancer[4] = 'mixed - ERROR'
                    del reps[minCheck[0]][0]

            elif len(minCheck)>1:
                if (initializedEnhancer == False):
                    thisEnhancer = reps[minCheck[0]][0]
                    initializedEnhancer = True

                if (int(reps[minCheck[0]][0][1])-int(thisEnhancer[2])>separation):
                    incompleteEnhancer = False
                    result.append(thisEnhancer)
                else:
                    for i in range(0,len(minCheck)):
                        if int(reps[minCheck[i]][0][2])>int(thisEnhancer[2]):
                            thisEnhancer[2]=reps[minCheck[i]][0][2]
                        if thisEnhancer[3] != reps[minCheck[i]][0][3]:
                            thisEnhancer[3] = 'mixed - ERROR'
                        if thisEnhancer[4] != reps[minCheck[i]][0][4]:
                            thisEnhancer[4] = 'mixed - ERROR'
                        del reps[minCheck[i]][0]

    return result

def loadFiles():
    fList = []
    process = False
    fileCount = 1
    while (process != True):
        try:
            f= open(sys.argv[fileCount],'r')
            fList.append(f)
            fileCount = fileCount + 1
        except:
            if fileCount==1:
                print('Must provide enhancers file(s) as argument.')
                print('Usage: python merge_enhancers.py [EHfile]*n')
                return process
            else:
                process = True
                print('Merging enhancers...')

    if (process == True):
        reps = [[]]*len(fList)

        #If + strand, filter out even motifs; if - strand, filter out odd motifs:
        #Remove "_random" chromosomes
        for i in range(0,len(fList)):
            reps[i]=[]
            current = fList[i].readline()
            current = current.strip().split('\t')
            while (current != ['']):
                if 'chrM' not in current[0]:
                    thisChrome = current[0]
                    currentChrome = []
                    while (current[0]==thisChrome and current != ['']):
                        currentChrome.append(current)
                        current = fList[i].readline()
                        current = current.strip().split('\t')
                    reps[i].append(currentChrome)
                else:
                    current = fList[i].readline()
                    current = current.strip().split('\t')
        return reps
    else:
        return False

def mergeJunctions():
    reps = loadFiles()

    largest = -float('inf')
    for i in range(0,len(reps)):
        if len(reps[i])>largest:
            largest = len(reps[i])

    for i in range(0,len(reps)):
        dif = largest - len(reps[i])
        for j in range(0,dif):
            reps[i].append([])

    if reps != False:
        merged = ['']*len(reps[0])
        for i in range(0,len(reps[0])):
            total = []
            for j in range(0,len(reps)):
                total.append(reps[j][i])
            merged[i]=kmerge(total)

        outputFile = open('enhancers_merged.bed','w+')
        for i in range(0, len(merged)):
            for j in range(0, len(merged[i])):
                k = '\t'.join(merged[i][j])
                outputFile.write(k +'\n')

    return True

if __name__ == "__main__":
#    t0 = time.time()
    merged = mergeJunctions()
#    t1 = time.time()
#    t2 = t1 - t0

#    if merged:
#        print 'Merged junctions: ' +str(t2) + ' seconds elapsed'
