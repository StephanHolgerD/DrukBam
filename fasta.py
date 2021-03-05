import pysam
import sys
import matplotlib.pyplot as plt

colorDict={'A':'red','C':'blue','T':'green','G':'yellow','-':'black','N':'pink'}
fig,ax=plt.subplots()
with pysam.FastaFile(sys.argv[1]) as fa:
    for _,n in enumerate(fa.fetch('EGFP', 10,20)):
        print(_)
        print(n)
        ax.text(10+_,1,n,color=colorDict[n])
ax.set(xlim=(10,20),ylim=(1,5))
plt.savefig('lol.pdf')
