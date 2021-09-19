import pysam
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from multiprocessing import Pool

class CalcMapping():
    def __init__(self,mapping,chrom,start,end,threads=1,coverage=200,flag='None',chunksize=1000):
        self.mapping=mapping
        self.chrom=chrom
        self.start=start-1
        self.end=end
        self.maxHeight=coverage
        self.Fontsize=3
        self.threads=threads
        self.flag = flag
        if self.flag is not None:
            print('flag')
            self.flag=int(self.flag)
            self.flags=set(self.SolveFlag(self.flag))
        self.chunksize=chunksize

    def SolveFlag(self,flag):
        allFlags=[1,2,4,8,16,32,64,128,256,512,1048]
        def OneLess(remainer,pot,FoundFlags=None):
            if FoundFlags is None:
                FoundFlags=[]
            if remainer ==0:
                return FoundFlags
            potF=[]
            for f in pot:
                if f == remainer:
                    FoundFlags.append(f)
                    return FoundFlags
                if f<remainer:
                    potF.append(f)
            remainer=remainer-max(potF)
            FoundFlags.append(max(potF))
            potF.remove(max(potF))
            if remainer == 0:
                return FoundFlags
            else:
                FoundFlags=OneLess(remainer,potF,FoundFlags)
                return FoundFlags
        return OneLess(flag,allFlags)

    def plotList(self,chrom,start,end,forwardOnly,reverseOnly):
        df=pd.DataFrame(0,index=range(0,self.maxHeight+50),columns=range(start,end+1))
        plotList=[]
        with pysam.AlignmentFile(self.mapping) as s:
            for record in tqdm(s.fetch(str(chrom),start,end,until_eof=True)):
            #    if self.flag is not None:
            #        print('lol')
            #        readFlags=set(self.SolveFlag(record.flag))
            #        if len(readFlags.intersection(self.flags))!=len(readFlags):
            #            continue
                #print(df)
                #print(list(df))
                if not record.is_unmapped:
                    if record.is_reverse:
                        if forwardOnly:
                            continue
                    if not record.is_reverse:
                        if reverseOnly:
                            continue
                
                #    if record.reference_start>=start and record.reference_end<end:
                #        if df.loc[self.maxHeight,start:end].sum()>=((self.end-self.start)*0.9):
                #            print('hardbreak')
                #            break

                #        if df.loc[self.maxHeight+50-1,record.reference_start:record.reference_end].sum()>=((record.reference_end-record.reference_start)*0.9):
                #            continue
                    endPos=record.reference_end
                    if record.reference_end>=end:
                        endPos=end
                    startPos=record.reference_start
                    if record.reference_start<=start:
                        startPos=start
                    if sum(df[range(startPos,endPos+1)].sum(axis=1)) == 0:
                        for _ in range(startPos,endPos+1):
                            if _ >=(end):
                                continue
                            df.at[0,_]=1
                        if record.is_reverse:

                            plotList.append((0,record.reference_start,record.reference_end,'r',
                                             record.qname+'_R',
                                             record.query_alignment_sequence,
                                             record.cigarstring,
                                             record.mate_is_unmapped))
                        if not record.is_reverse:

                            plotList.append((0,record.reference_start,record.reference_end,'f',
                                             record.qname+'_F',
                                             record.query_alignment_sequence,
                                             record.cigarstring,
                                             record.mate_is_unmapped))


                    else:

                        for p,v in enumerate(list(df[range(startPos,endPos+1)].sum(axis=1))):
                            if v ==0:
                                for _ in range(startPos,endPos+1):
                                    df.at[p,_]=1

                                if record.is_reverse:

                                    plotList.append((p,record.reference_start,record.reference_end,'r',
                                                     record.qname+'_R',
                                                     record.query_alignment_sequence,
                                                     record.cigarstring,
                                                     record.mate_is_unmapped))
                                if not record.is_reverse:

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
                    for _ in range(record.reference_start,record.reference_start+self.chunk_size+100):
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
