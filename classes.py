import pysam
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from multiprocessing import Pool
        
class PlotMapping():
    def __init__(self,mapping,chrom,start,end):
        self.mapping=mapping
        self.chrom=chrom
        self.start=start-1
        self.end=end
        self.maxHeight=100
        self.Fontsize=3
        
    
        
    def plotList(self,chrom,start,end):
        df=pd.DataFrame(0,index=range(0,1000),columns=range(start,end+1))
        plotList=[]
        with pysam.AlignmentFile(self.mapping) as s:
            for record in tqdm(s.fetch(str(chrom),start,end,until_eof=True)):
                if not record.is_unmapped:
                    for _ in range(record.reference_start,record.reference_end+1):
                        if _ not in list(df):
                            df[_]=0
                    if sum(df[range(record.reference_start,record.reference_end+1)].sum(axis=1)) == 0:
                        for _ in range(record.reference_start,record.reference_end+1):
                            df.at[0,_]=1
                        if record.is_reverse:
                            plotList.append((0,record.reference_start,record.reference_end,'r',
                                             record.qname+'_R',record.query_alignment_sequence,record.cigarstring))
                        else:
                            plotList.append((0,record.reference_start,record.reference_end,'f',
                                             record.qname+'_F',record.query_alignment_sequence,record.cigarstring))
                            
                    
                    else:
                        for p,v in enumerate(list(df[range(record.reference_start,record.reference_end+1)].sum(axis=1))):
                            if v ==0:
                                for _ in range(record.reference_start,record.reference_end+1):
                                    df.at[p,_]=1
                                if record.is_reverse:
                                    plotList.append((p,record.reference_start,record.reference_end,'r',
                                                     record.qname+'_R',record.query_alignment_sequence,record.cigarstring))
                                else:
                                    plotList.append((p,record.reference_start,record.reference_end,'f',
                                                     record.qname+'_F',record.query_alignment_sequence,record.cigarstring))
                                break
        return plotList
    
    def PlotreadsDF(self,lbefore,l,end_prevChunk):
        before=pd.DataFrame(lbefore,columns=['y','start','end','direction','name','qSeq','cigar'])
        before=before[before['end']>end_prevChunk]
        recent=pd.DataFrame(l,columns=['y','start','end','direction','name','qSeq','cigar'])
        for overlappRead in set(before['name']).intersection(set(recent['name'])):
            ybefore=before[before['name']==overlappRead]['y'].values[0]
            beforeReads=before[before['y']==ybefore]['name']
            yrecent=recent[recent['name']==overlappRead]['y'].values[0]
            recentReads=recent[recent['y']==yrecent]['name']
            recenUpdateReads=recent[recent['y']==ybefore]['name']
        
            for n in recenUpdateReads:
                i=recent[recent['name']==n].index.values[0]
                recent.at[i,'y']=yrecent
        
            for n in recentReads:
                i=recent[recent['name']==n].index.values[0]
                recent.at[i,'y']=ybefore
        return recent
    def startPlot(self,cols):
        if cols==1:
            print(self.end-self.start)
            if self.end-self.start <=100:
                x=0.07*(self.end-self.start)
                y=0.12*self.maxHeight
                fig,ax=plt.subplots(1,cols,figsize=(x,y))
                self.Fontsize=6
            else:
                print('jo')
                x=0.035*(self.end-self.start)
                print(x)
                y=0.06*self.maxHeight
                fig,ax=plt.subplots(1,cols,figsize=(x,y))
        if cols >1:
            x=0.035*5000
            y=0.06*self.maxHeight
            fig,ax=plt.subplots(1,cols,figsize=(x*cols,6))
            
        self.ax=ax
        self.fig=fig
        plt.subplots_adjust(wspace=0)
        
        
    def plotChunk(self,df,ax,start,end):
        for y,s,e,d in zip(df['y'],df['start'],df['end'],df['direction']):
            e=e
            s=s+1
            if d=='r':
                ax.plot((s,e),(y,y),color='grey')
                ax.plot((s+1,e-1),(y,y),linewidth=0.1,color='white')
                
            else:
                ax.plot((s,e),(y,y),color='black')
                ax.plot((s+1,e-1),(y,y),linewidth=0.1,color='white')
                
        ax.set(xlim=(start,end),ylim=(0,100))
        ax.get_yaxis().set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
    
    def plotEmptyChunk(self,ax,start,end):
        ax.set(xlim=(start,end),ylim=(0,100))
        ax.get_yaxis().set_visible(False)
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
    
    
    def PlotNucChunk(self,df,ax,start,end):
        colorDict={'A':'red','C':'blue','T':'green','G':'yellow','-':'black','N':'pink'}
        for y,s,e,d,qS,cig in zip(df['y'],df['start'],df['end'],df['direction'],df['qSeq'],df['cigar']):
            e=e
            s=s+1
            if y>100:
                continue
            chunk_cigarstring=self.CigChunker(cig)
            query_alignment_sequence=qS
            for p,_ in enumerate(chunk_cigarstring):
                if _=='D':
                    qs1=query_alignment_sequence[:p]
                    qs2=query_alignment_sequence[p:]
                    query_alignment_sequence=qs1+'-'+qs2
                    
                    
            chunk_cigarstring=[x for x in chunk_cigarstring if x !='S']
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
                            color='black',family='monospace',ha='center',va='center',
                            bbox=dict(boxstyle='square,pad=0', fc=colorDict[query_alignment_sequence[p]], ec='none'))
                    continue
                if alignPos=='I':
                    ax.plot((x-1,x-1),(y-0.5,y+0.5),linewidth=0.5)
                    continue
                if alignPos=='D':
                    ax.text(x,y,query_alignment_sequence[p], fontsize=self.Fontsize,ha='center',va='center',
                            color=colorDict[query_alignment_sequence[p]])
                    continue
    
    def Plot(self):
        span=self.end-self.start
        multi=[]
        start=self.start
        end=self.start+1000
        for _ in range(0,int(span/1000)):
            multi.append((self.chrom,start,end))
            start=end
            end=end+1000
        if span/1000 - int(span/1000)!=0:
            multi.append((self.chrom,start,int(start+((span/1000 - int(span/1000))*1000))))
        with Pool(processes=3) as pool:
            results = pool.starmap(self.plotList, multi)
        self.startPlot(len(results))
        for p,chunk in enumerate(zip(results,multi)):
            
            if results[p]==[] and len(results)==1:
                self.plotEmptyChunk(self.ax,chunk[1][1],chunk[1][2])
                continue
            
            if results[p]==[]:
                self.plotEmptyChunk(self.ax[p],chunk[1][1],chunk[1][2])
            
            
            if len(results)==1:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar'])
                self.plotChunk(d,self.ax,chunk[1][1],chunk[1][2])
                self.PlotNucChunk(d,self.ax,chunk[1][1],chunk[1][2])
                continue
                
            if p==0:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar'])
                self.plotChunk(d,self.ax[p],chunk[1][1],chunk[1][2])
            else:
                if results[p-1]==[]:
                    d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar'])
                    self.plotChunk(d,self.ax[p],chunk[1][1],chunk[1][2])
                else:
                    d=self.PlotreadsDF(results[p-1],results[p],multi[p-1][2])
                    self.plotChunk(d,self.ax[p],chunk[1][1],chunk[1][2])
        plt.savefig('{}_{}_{}.pdf'.format(self.chrom,str(self.start),str(self.end)))
                    

