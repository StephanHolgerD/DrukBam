import pysam
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from multiprocessing import Pool

class CalcPlot():
    def __init__(self,mapping,chrom,start,end,threads=1,coverage=200,flag='None',chunksize=1000):
        self.mapping=mapping
        self.chrom=chrom
        self.start=start-1
        self.end=end
        self.maxHeight=coverage
        self.Fontsize=3
        self.threads=threads
        self.flag=flag
        self.chunksize=chunksize

    def startPlot(self,cols,direction,schematic):

        self.cols=cols
        if not direction and not schematic:
            if cols==1:
                if self.end-self.start <=100:
                    x=0.07*(self.end-self.start)
                    y=0.12*self.maxHeight
                    fig,ax=plt.subplots(1,cols,figsize=(x,y))
                    self.Fontsize=6
                else:
                    x=0.035*(self.end-self.start)
                    y=0.06*self.maxHeight
                    fig,ax=plt.subplots(1,cols,figsize=(x,y))
            if cols >1:
                x=0.035*5000
                y=0.06*self.maxHeight
                fig,ax=plt.subplots(1,cols,figsize=(x*cols,6))
        if direction and not schematic:
            if cols==1:
                if self.end-self.start <=100:
                    x=0.07*(self.end-self.start)
                    y=0.12*self.maxHeight*2
                    fig,ax=plt.subplots(2,cols,figsize=(x,y))
                    self.Fontsize=6
                else:
                    x=0.035*(self.end-self.start)
                    y=0.06*self.maxHeight*2
                    fig,ax=plt.subplots(2,cols,figsize=(x,y))
            if cols >1:
                x=0.035*5000
                y=0.06*self.maxHeight*2
                fig,ax=plt.subplots(2,cols,figsize=(x*cols,6))


        if not direction and schematic:
            if self.cols >= 5:
                y=(0.12*self.maxHeight)/5
                x=12
                self.innerLW=0
                self.LW=0.2
            else:
                self.innerLW=0.4/self.cols
                self.LW=3/self.cols
                x=12
                y=(0.12*self.maxHeight)/self.cols
            fig,ax=plt.subplots(1,self.cols,figsize=(x,y))
        if  direction and schematic:
            if self.cols >= 5:
                y=((0.12*self.maxHeight)/5)*2
                x=12
                self.innerLW=0
                self.LW=0.2
            else:
                self.innerLW=0.4/self.cols
                self.LW=3/self.cols
                x=12
                y=((0.12*self.maxHeight)/self.cols)*2
            fig,ax=plt.subplots(2,self.cols,figsize=(x,y))

        plt.suptitle('chromosome {} from {} to {}'.format(str(self.chrom),str(self.start),str(self.end)))
        self.ax=ax
        self.fig=fig
        plt.subplots_adjust(wspace=0,hspace=0.05)
        plt.style.use('bmh')


    def plotChunk(self,df,ax,start,end,direction=None):
        df['y']=df['y']+1

        for y,s,e,d in zip(df['y'],df['start'],df['end'],df['direction']):
            e=e
            s=s+1
            if y>(self.maxHeight):
                ax.plot((s,e),(self.maxHeight+1,self.maxHeight+1),color='red',alpha=0.1)
                continue
            if d=='r':
                ax.plot((s,e),(y,y),color='grey',linewidth=self.LW)
                ax.plot((s+1,e-1),(y,y),linewidth=self.innerLW,color='white')

            else:
                ax.plot((s,e),(y,y),color='black',linewidth=self.LW)
                ax.plot((s+1,e-1),(y,y),linewidth=self.innerLW,color='white')



    def AxSet(self,ax,start,end,chunk=False,direction='all'):
        ylabel=direction
        if not chunk:
            ax.set(xlim=(start,end),ylim=(0,self.maxHeight+2),ylabel=ylabel)
            ax.set_xticks([start,end],)
            ax.set_xticklabels([str(start/1_000_000)+' mb',str(end/1_000_000)+' mb'], rotation=40, ha='right')

        else:
            ax.get_yaxis().set_visible(False)
            ax.set(xlim=(start,start+self.chunksize),ylim=(0,self.maxHeight+2))
            if start==self.start:
                ax.set(ylabel=ylabel)
                ax.get_yaxis().set_visible(True)
                ax.set_xticks([start])
                ax.set_xticklabels([str(start/1_000_000)+' mb'], rotation=40, ha='right',fontsize=6)
            if end==self.end:
                ax.set_xticks([start,start+self.chunksize])
                ax.set_xticklabels(['',str((start+self.chunksize)/1_000_000)+' mb'], rotation=40, ha='right',fontsize=6)
            if  start!=self.start and end!=self.end:
                ax.set_xticks([start])
                ax.set_xticklabels([''], rotation=40, ha='right',fontsize=6)


        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)

    def CigChunker(self,cig):
        cigL=[]
        def parseCig(cig):
            for p,l in enumerate(cig):
                if l.isalpha():
                    liste=([cig[p]]*int(cig[:p]))
                    cig=cig[p+1:]
                    return cig,liste
        while cig:
            cigL=cigL+parseCig(cig)[1]
            cig=parseCig(cig)[0]
        return cigL


    def PlotNucChunk(self,df,ax,start,end,flag='None'):
        df['y']=df['y']+1
        colorDict={'A':'red','C':'blue','T':'green','G':'yellow','-':'black','N':'pink'}
        for y,s,e,d,qS,cig,mate in zip(df['y'],df['start'],df['end'],df['direction'],df['qSeq'],df['cigar'],df['mateMap']):
            e=e
            s=s+1
            if y>self.maxHeight:
                ax.plot((s,e),(self.maxHeight+1,self.maxHeight+1),color='red',alpha=0.1)
                continue
            chunk_cigarstring=self.CigChunker(cig)
            query_alignment_sequence=qS
            for p,_ in enumerate(chunk_cigarstring):
                if _=='D':
                    qs1=query_alignment_sequence[:p]
                    qs2=query_alignment_sequence[p:]
                    query_alignment_sequence=qs1+'-'+qs2

            chunkL=len(chunk_cigarstring)
            chunk_cigarstringS=[x for x in chunk_cigarstring if x =='S']
            chunk_cigarstring=[x for x in chunk_cigarstring if x !='S']
            softClipp=len(chunk_cigarstringS)/chunkL
            alpha=1
            if flag=='None':
                alpha=1
            if flag=='MateUnmapped':
                if mate:
                    alpha=0.3
            if flag=='SoftClipped':
                if softClipp>=0.05:
                    alpha=0.3
            if flag=='MateUnmappedSoftClipped':
                if softClipp>=0.1 or mate:
                    alpha=0.3
            for p,alignPos in enumerate(chunk_cigarstring):
                x=s+p
                if x>self.end:
                    continue
                if x<self.start:
                    continue
                if alignPos=='S':
                    continue
                if alignPos=='M':
                    ax.text(x,y,query_alignment_sequence[p],fontsize=self.Fontsize,
                            color='black',alpha=alpha,family='monospace',ha='center',va='center',
                            bbox=dict(alpha=alpha,boxstyle='square,pad=0', fc=colorDict[query_alignment_sequence[p]], ec='none'))
                    continue
                if alignPos=='I':
                    ax.plot((x-1,x-1),(y-0.5,y+0.5),linewidth=0.5)
                    continue
                if alignPos=='D':
                    ax.text(x,y,query_alignment_sequence[p], fontsize=self.Fontsize,ha='center',va='center',
                            color=colorDict[query_alignment_sequence[p]])
                    continue
