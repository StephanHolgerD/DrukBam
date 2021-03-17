import pysam
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from multiprocessing import Pool

class CalcMapping():
    def __init__(self,mapping,chrom,start,end,threads=1,coverage=200,flag='None'):
        self.mapping=mapping
        self.chrom=chrom
        self.start=start-1
        self.end=end
        self.maxHeight=coverage
        self.Fontsize=3
        self.threads=threads
        self.flag=flag

    def plotList(self,chrom,start,end):
        df=pd.DataFrame(0,index=range(0,self.maxHeight+50),columns=range(start,end+1))
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
                                             record.qname+'_R',
                                             record.query_alignment_sequence,
                                             record.cigarstring,
                                             record.mate_is_unmapped))
                        else:
                            plotList.append((0,record.reference_start,record.reference_end,'f',
                                             record.qname+'_F',
                                             record.query_alignment_sequence,
                                             record.cigarstring,
                                             record.mate_is_unmapped))


                    else:
                        for p,v in enumerate(list(df[range(record.reference_start,record.reference_end+1)].sum(axis=1))):
                            if v ==0:
                                for _ in range(record.reference_start,record.reference_end+1):
                                    df.at[p,_]=1

                                if record.is_reverse:
                                    plotList.append((p,record.reference_start,record.reference_end,'r',
                                                     record.qname+'_R',
                                                     record.query_alignment_sequence,
                                                     record.cigarstring,
                                                     record.mate_is_unmapped))
                                else:
                                    plotList.append((p,record.reference_start,record.reference_end,'f',
                                                     record.qname+'_F',
                                                     record.query_alignment_sequence,
                                                     record.cigarstring,
                                                     record.mate_is_unmapped))
                                break
        return plotList

    def plotListRVRS(self,chrom,start,end):
        df=pd.DataFrame(0,index=range(0,self.maxHeight+50),columns=range(start,end+1))
        plotList=[]
        with pysam.AlignmentFile(self.mapping) as s:
            for record in tqdm(s.fetch(str(chrom),start,end,until_eof=True)):
                if not record.is_unmapped and record.is_reverse:
                    #print('s')
                    #print(record.reference_start)
                    #record.reference_start=record.reference_start+1
                    #record.reference_end=record.reference_end+1
                    #print(record.reference_start)
                    #print('e')
                    for _ in range(record.reference_start,record.reference_end):
                        if _ not in list(df):
                            df[_]=0
                    if sum(df[range(record.reference_start,record.reference_end)].sum(axis=1)) == 0:
                        for _ in range(record.reference_start,record.reference_end):
                            df.at[0,_]=1

                        plotList.append((0,record.reference_start,record.reference_end,'r',
                                             record.qname+'_R',
                                             record.query_alignment_sequence,
                                             record.cigarstring,
                                             record.mate_is_unmapped))


                    else:
                        for p,v in enumerate(list(df[range(record.reference_start,record.reference_end)].sum(axis=1))):
                            if v ==0:
                                for _ in range(record.reference_start,record.reference_end):
                                    df.at[p,_]=1
                                plotList.append((p,record.reference_start,record.reference_end,'r',
                                                     record.qname+'_R',
                                                     record.query_alignment_sequence,
                                                     record.cigarstring,
                                                     record.mate_is_unmapped))
                                break
        return plotList

    def plotListFRWRD(self,chrom,start,end):
        df=pd.DataFrame(0,index=range(0,self.maxHeight+50),columns=range(start,end+1))
        plotList=[]
        with pysam.AlignmentFile(self.mapping) as s:
            for record in tqdm(s.fetch(str(chrom),start,end,until_eof=True)):
                if not record.is_unmapped and not record.is_reverse:
                    for _ in range(record.reference_start,record.reference_end+1):
                        if _ not in list(df):
                            df[_]=0
                    if sum(df[range(record.reference_start,record.reference_end+1)].sum(axis=1)) == 0:
                        for _ in range(record.reference_start,record.reference_end+1):
                            df.at[0,_]=1

                        plotList.append((0,record.reference_start,record.reference_end,'f',
                                             record.qname+'_F',
                                             record.query_alignment_sequence,
                                             record.cigarstring,
                                             record.mate_is_unmapped))


                    else:
                        for p,v in enumerate(list(df[range(record.reference_start,record.reference_end+1)].sum(axis=1))):
                            if v ==0:
                                for _ in range(record.reference_start,record.reference_end+1):
                                    df.at[p,_]=1
                                plotList.append((p,record.reference_start,record.reference_end,'f',
                                                     record.qname+'_F',
                                                     record.query_alignment_sequence,
                                                     record.cigarstring,
                                                     record.mate_is_unmapped))
                                break
        return plotList



    def PlotreadsDF(self,lbefore,l,end_prevChunk):
        before=pd.DataFrame(lbefore,columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
        before=before[before['end']>end_prevChunk]
        recent=pd.DataFrame(l,columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
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
