import sys
import os
import pandas as pd
import numpy as np
import matplotlib
if os.environ.get('DISPLAY','') == '':
    print('No display found. Using non-interactive Agg backend')
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns



tag = sys.argv[1]

enhfl = tag+'_arranged.txt'

cols = ['name','nm','chrom','strand','start','end']
df = pd.read_csv(enhfl,sep='\t',header=None,names = cols)
df=df.set_index('nm')

enhancer_lengths=np.asarray(list(df['end']-df['start']))
avg_length=np.mean(enhancer_lengths)
size_infig=round(avg_length,2)

sns_plot=sns.distplot(enhancer_lengths).set_title('Enhancer Size Distribution')
plt.figtext(.6, .6, 'mean = '+str(size_infig))
plt.yticks([])
plt.show()
fig = sns_plot.get_figure()
fig.savefig(tag+"_dist-size.png")

dropfl = open(tag+'_dropped.txt','w')

tlist = df.values.tolist()
insigs = []
filtered_df=df
count_lows = 0
for i in range(0,len(enhancer_lengths)):
    if enhancer_lengths[i]<=25:
        count_lows = count_lows + 1
        insigs.append(i+1)
        filtered_df=filtered_df.drop('ENH_NM_'+str(i+1))
        out='\t'.join(map(str,tlist[i]))
        dropfl.write(out+'\tinsignificant length (<=25bp)\n')

dropfl.close()

filtered_df=filtered_df.reset_index()
bed_cols = ['chrom','start','end','name','nm','strand']
new_txt=filtered_df[cols]
new_bed=filtered_df[bed_cols]

new_txt.to_csv(tag+'_glist.txt',sep='\t',index=False)
new_bed.to_csv(tag+'_glist.bed',sep='\t',index=False,header=False)
