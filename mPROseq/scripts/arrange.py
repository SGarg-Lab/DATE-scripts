import sys

c = 1
with open(sys.argv[1]) as f:
    for inline in f:
        ln = inline.strip().split('\t')
        #score = round(float(ln[3]),3)  ###score derived from dREG, not dREG.HD
        #print('enhancer.'+str(c)+'_'+str(score)+'\t'+'ENH_NM_'+str(c)+'\t'+ln[0]+'\t+\t'+ln[1]+'\t'+ln[2])
        print('enhancer.'+str(c)+'\t'+'ENH_NM_'+str(c)+'\t'+ln[0]+'\t+\t'+ln[1]+'\t'+ln[2])
        c = c + 1
