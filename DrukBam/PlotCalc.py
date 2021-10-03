import pysam
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from multiprocessing import Pool
import configparser
import os
import sys
import matplotlib.patches as patches
class CalcPlot():
    def __init__(self,mapping,chrom,start,end,threads=1,coverage=200,flag='None',chunksize=1000, fasta=None,style='classic',outlineoff=False,clipped=1):
        self.mapping=mapping
        self.chrom=chrom
        self.start=start
        self.end=end
        self.maxHeight=coverage
        self.Fontsize=3
        self.threads=threads
        self.flag=flag
        self.chunksize=chunksize
        self.fasta = fasta
        self.outlineoff=outlineoff
        config = configparser.ConfigParser()
        if style!=None:
            if os.path.isfile(style):
                config.read(style)
            else:
                print('{} is not a existing file'.format(style))
                config.read(os.path.join(os.path.dirname(__file__),'classic.ini'))

        else:
            config.read(os.path.join(os.path.dirname(__file__),'classic.ini'))
        self.clipped=clipped
        colorDict={}

        for a in config['nucleotide color']:
            colorDict[a.upper()]=config['nucleotide color'][a]
            colorDict[a.lower()]=config['nucleotide color'][a]
        for a in config['special chars']:
            colorDict[a]=config['special chars'][a]
        for a in config['MatplotStyle']:
            colorDict[a]=config['MatplotStyle'][a]
        colorDict['ClippedAlpha']=float(config['clipped']['alpha'])
        colorDict['ClippedBackground']=config['clipped']['background']
        colorDict['ClippedColor']=config['clipped']['outline']

        colorDict['OutlineAlpha']=float(config['outline']['alpha'])
        colorDict['OutlineBackground']=config['outline']['background']
        colorDict['OutlineColor']=config['outline']['outline']

        colorDict['xticks']=config['Figure']['tickspacing']

        self.colorDict=colorDict
        plt.style.use(self.colorDict['pltstyle'])

        plt.rc('font', size=8)          # controls default text sizes
        plt.rc('axes', titlesize=8)     # fontsize of the axes title
        plt.rc('axes', labelsize=8)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=6)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=6)    # fontsize of the tick labels
        plt.rc('legend', fontsize=8)    # legend fontsize
        plt.rc('figure', titlesize=8)  # fontsize of the figure title



        #self.colorDict={'A':'red','C':'blue','T':'green','G':'yellow','-':'white','N':'pink','a':'red','c':'blue','t':'green','g':'yellow','-':'white','n':'pink','|':'black'}
    def startPlot(self,cols,direction,schematic):
        if self.end-self.start <=100:
            self.Fontsize=self.Fontsize/2
            x_perN=0.07/2
            y_perN=0.12/2
            self.LW=0.1
            self.Fontsize=3
            self.OutlineInnerW=2.2
            self.OutlineOuterW=2.6



        else:
            x_perN=0.07
            y_perN=0.12
            self.LW=0.2
            self.Fontsize=3*2
            self.OutlineInnerW=2.2*2
            self.OutlineOuterW=2.6*2

        msize=.2
        self.cols=cols
        self.innerLW=0
        if not direction and not schematic:
            if cols==1:
                x=x_perN*(self.end-self.start)
                y=y_perN*self.maxHeight
                self.markersize=msize
                fig,ax=plt.subplots(1,cols,figsize=(x,y))
            if cols >1:
                x=x_perN*5000
                y=y_perN*self.maxHeight
                fig,ax=plt.subplots(1,cols,figsize=(x*cols,6))
        if direction and not schematic:
            if cols==1:
                if self.end-self.start <=100:
                    x=x_perN*(self.end-self.start)
                    y=y_perN*self.maxHeight*2
                    self.markersize=msize
                    fig,ax=plt.subplots(2,cols,figsize=(x,y))
                else:
                    x=x_perN*(self.end-self.start)
                    y=y_perN*self.maxHeight*2
                    self.markersize=msize
                    fig,ax=plt.subplots(2,cols,figsize=(x,y))
            if cols >1:
                x=x_perN*5000
                y=y_perN*self.maxHeight*2
                fig,ax=plt.subplots(2,cols,figsize=(x*cols,6))


        if not direction and schematic:
            if self.cols >= 5:
                y=(y_perN*self.maxHeight)/5
                x=12
                self.innerLW=0
                self.LW=0.2
            else:
                self.innerLW=0.4/self.cols
                self.LW=2/self.cols
                x=12
                y=(y_perN*self.maxHeight)/self.cols
            fig,ax=plt.subplots(1,self.cols,figsize=(x,y))
        if  direction and schematic:
            if self.cols >= 5:
                y=((y_perN*self.maxHeight)/5)*2
                x=12
                self.innerLW=0
                self.LW=0.2
            else:
                self.innerLW=0.4/self.cols
                self.LW=2/self.cols
                x=12
                y=((y_perN*self.maxHeight)/self.cols)*2
            fig,ax=plt.subplots(2,self.cols,figsize=(x,y))

        plt.suptitle('{} {}-{}'.format(str(self.chrom),str("{:,}".format(self.start)),str("{:,}".format(self.end))))
        self.ax=ax
        self.fig=fig
        plt.subplots_adjust(wspace=0,hspace=0.05, )
        return self.fig,self.ax


    def plotChunk(self,df,ax,start,end,direction=None):
        df['y']=df['y']+1

        for y,s,e,d in zip(df['y'],df['start'],df['end'],df['direction']):
            e=e
            s=s
            if y>(self.maxHeight):
                ax.plot((s,e),(self.maxHeight+1,self.maxHeight+1),color='red',alpha=0.1)
                continue
            if d=='r':
                ax.plot((s,e),(y,y),color='grey',linewidth=self.LW)
                ax.plot((s+1,e-1),(y,y),linewidth=self.innerLW,color='white')

            else:
                ax.plot((s,e),(y,y),color='black',linewidth=self.LW)
                ax.plot((s+1,e-1),(y,y),linewidth=self.innerLW,color='white')





    def Xspacing(self,start,end,size=100,vcfPos=False):
        if int(end-start)<size:
            if not type(vcfPos)==int:
                return [start,end]
            else:
                return [start,vcfPos,end]

        digitDict={100:2,1000:3,10:1,10_000:4}
        diffStart=size-int(str(start)[-1*digitDict[size]:])
        rangeStart=start+diffStart
        diffEnd=int(str(end)[-1*digitDict[size]:])
        rangeEnd=end-diffEnd
        x=size
        if rangeEnd==end:
            x=0
        sizeChunks=list(range(rangeStart,rangeEnd+x,size))
        chunks=[start]+sizeChunks+[end]
        if type(vcfPos)==int:
            print('here')
            chunks=chunks+[vcfPos]
            chunks.sort()
        return chunks



    def AxSet(self,ax,start,end,chunk=False,direction='all',vcf=False):

        ylabel=direction
        if not chunk:
            xlim=self.Xspacing(start,end,size=int(self.colorDict['xticks']))
            ax.set(xlim=(start,end),ylim=(0,self.maxHeight+2),ylabel=ylabel)
            ax.yaxis.set_ticks_position('left')

            if vcf:
                xlim=self.Xspacing(start,end,size=int(self.colorDict['xticks']),vcfPos=int(start+(end-start)/2))
                ax.set_xticks(xlim)
                xlimStr=[str("{:,}".format(x)) for x in xlim]
                ax.set_xticklabels(xlimStr, rotation=40, ha='right')
            else:
                xlim=self.Xspacing(start,end,size=int(self.colorDict['xticks']))
                ax.set_xticks(xlim)
                xlimStr=[str("{:,}".format(x)) for x in xlim]
                ax.set_xticklabels(xlimStr, rotation=40, ha='right')

        else:
            ax.xaxis.set_tick_params(length=1)
            ax.get_yaxis().set_visible(False)
            xlim=self.Xspacing(start,start+self.chunksize,size=int(self.colorDict['xticks']))
            xlimStr=[str("{:,}".format(x)) for x in xlim]
            ax.set(xlim=(start,start+self.chunksize),ylim=(0,self.maxHeight+2))
            ax.yaxis.set_ticks_position('none')
            ax.xaxis.set_ticks_position('bottom')
            if start==self.start:
                ax.set(ylabel=ylabel)
                ax.get_yaxis().set_visible(True)
                ax.set_xticks(xlim[:-1])
                ax.set_xticklabels(xlimStr[:-1], rotation=40, ha='right')
            if end==self.end:
                ax.set_xticks(xlim[1:])
                ax.set_xticklabels(xlimStr[1:], rotation=40, ha='right')
            if  start!=self.start and end!=self.end:
                ax.set_xticks(xlim[1:-1])
                ax.set_xticklabels(xlimStr[1:-1], rotation=40, ha='right')


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



    def PlotFasta(self,ax):
        if self.fasta=='None':
            return
        with pysam.FastaFile(self.fasta) as fa:
            for _,n in enumerate(fa.fetch(self.chrom,self.start-1,self.end)):
                ax.text(self.start+_,0,n,fontsize=self.Fontsize,
                        color=self.colorDict[n],alpha=1,family='monospace',ha='center',va='center',
                        bbox=dict(alpha=1,boxstyle='square,pad=0', fc=self.colorDict['fasta background'], ec='none'))

    def listConsec(self,l):
        retL=[]
        c=0
        for e,x in enumerate(l):
            if x+1 not in l:
                retL.append(l[c:e+1])
                c=e+1
        return(retL)
    def PlotNucChunk(self,df,ax,start,end,flag='None'):
        plotC=set()
        #if self.readoutline:
    #        self.plotChunk(df,ax,self.start,self.end)
        df['y']=df['y']+1
        for y,s,e,d,qS,cig,mate in zip(df['y'],df['start'],df['end'],df['direction'],df['qSeq'],df['cigar'],df['mateMap']):
            if self.fasta != 'None':
                with pysam.FastaFile(self.fasta) as fa:
                    fastaChunk=str(fa.fetch(self.chrom,s,e)).upper()
            e=e+1
            s=s+1
            if y>self.maxHeight:
                ax.plot((s,e),(self.maxHeight+1,self.maxHeight+1),color=self.colorDict['max coverage'],alpha=0.1)
                continue

            chunk_cigarstring=self.CigChunker(cig)
            query_alignment_sequence=qS
            chunkL=len(chunk_cigarstring)
            chunk_cigarstringS=[x for x in chunk_cigarstring if x =='S' or x =='H']
            chunk_cigarstring=[x for x in chunk_cigarstring if x !='S' and x !='H']

            ipos=[x for x,y in enumerate(chunk_cigarstring) if y=='I']
            ipos=self.listConsec(ipos)
            chunk_cigarstring=[x for x in chunk_cigarstring if x !='I']
            iposCounter=0

            for i in ipos:
                query_alignment_sequence="".join([x for _,x in enumerate(query_alignment_sequence) if _+iposCounter not in i])
                iposCounter=iposCounter+len(i)

            for p,_ in enumerate(chunk_cigarstring):
                if _=='D' or _=='N':
                    qs1=query_alignment_sequence[:p]
                    qs2=query_alignment_sequence[p:]
                    query_alignment_sequence=qs1+'-'+qs2


            softClipp=len(chunk_cigarstringS)/chunkL

            alpha=self.colorDict['OutlineAlpha']
            outlineColor=self.colorDict['OutlineColor']
            outlineBackground=self.colorDict['OutlineBackground']

            if softClipp >= self.clipped:
                alpha=self.colorDict['ClippedAlpha']
                outlineColor=self.colorDict['ClippedColor']
                outlineBackground=self.colorDict['ClippedBackground']



            fastapos=0
            drawChunk=[]
            if not self.outlineoff:
                ax.plot((s,s+len(query_alignment_sequence)-1),(y,y),color=outlineColor,linewidth=self.OutlineOuterW,zorder=0)
                ax.plot((s+0.1,s+len(query_alignment_sequence)-1-0.1),(y,y),color=outlineBackground,linewidth=self.OutlineInnerW,alpha=alpha,zorder=0)

            fontDict={'A':r'$\mathtt{A}$','C':r'$\mathtt{C}$','G':r'$\mathtt{G}$','T':r'$\mathtt{T}$','-':r'$\mathtt{-}$'}
            xs=[x+s for x in range(0,len(chunk_cigarstring))]
            ys=len(chunk_cigarstring)*[y]
            if self.fasta=='None':
                for _ in ['A','C','T','G','-']:
                    xp=[x for x,nuc in zip(xs,query_alignment_sequence) if nuc == _]
                    yp=[y for y,nuc in zip(ys,query_alignment_sequence) if nuc == _]
                    ax.scatter(xp,yp,marker=fontDict[_],lw=0,color=self.colorDict[_],s=4,alpha=alpha)
                continue

            if self.fasta != 'None':
                ax.scatter([x+s for x in range(0,len(chunk_cigarstring))],len(chunk_cigarstring)*[y], marker='.',s=.2,linewidth=self.LW,color=self.colorDict['dot'],alpha=alpha)

            for p,alignPos in enumerate(chunk_cigarstring):
                x=s+p
                if x>self.end:
                    continue
                if x<self.start:
                    fastapos=fastapos+1
                    continue
                if alignPos=='S':
                    continue
                if alignPos=='N':
                    ax.text(x,y,query_alignment_sequence[p], fontsize=self.Fontsize,ha='center',va='center',family='monospace',
                            color=self.colorDict['-'],bbox=dict(alpha=alpha,boxstyle='square,pad=0', fc=self.colorDict['gap background'], ec='darkgrey',linewidth=0.00001))
                    fastapos=fastapos+1
                    continue
                if alignPos=='M':
                    if self.fasta != 'None' and fastaChunk[fastapos]==query_alignment_sequence[p]:
                        fastapos=fastapos+1
                        continue
                    else:
                        if self.fasta != 'None':
                            ax.text(x,y,query_alignment_sequence[p],fontsize=self.Fontsize,
                                color=self.colorDict['nuc missmatch font'],alpha=alpha,family='monospace',ha='center',va='center_baseline',zorder=2
                                )
                                #bbox=dict(alpha=alpha,boxstyle='square,pad=0.1' ,fc=self.colorDict[query_alignment_sequence[p]], ec='none')
                            rect = patches.Rectangle((x-0.45, y-0.5),
                            1,1,edgecolor='none',facecolor=self.colorDict[query_alignment_sequence[p]],zorder=1)
# add rectangle to plot
                            ax.add_patch(rect)
                            fastapos=fastapos+1
                            continue

                if alignPos=='D':
                    ax.text(x,y,query_alignment_sequence[p], fontsize=self.Fontsize,ha='center',va='center',family='monospace',
                            color=self.colorDict['-'],bbox=dict(alpha=alpha,boxstyle='square,pad=0', fc=self.colorDict['gap background'], ec='darkgrey',linewidth=0.00001))
                    fastapos=fastapos+1
                    continue

            if ipos!=[]:
                for i in ipos:
                    minimum=min(i)
                    x=s+minimum
                    for ii in i:
                        ax.plot((x,x),(y-0.5,y+0.5),linewidth=0.5, color=self.colorDict['insertion'])
                        x=x+0.05
