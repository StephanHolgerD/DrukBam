import pysam
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from multiprocessing import Pool
from classes.bamCalc import CalcMapping
from classes.PlotCalc import CalcPlot
import sys
class PlotMapping():
    def __init__(self,mapping,chrom,start,end,threads=1,coverage=200,flag='None',direction=False,schematic=False,chunksize=1000):
        self.mapping=mapping
        self.chrom=chrom
        self.start=start-1
        self.end=end
        self.maxHeight=coverage
        self.Fontsize=3
        self.threads=threads
        self.flag=flag
        self.chunksize=chunksize
        self.CalcMapping=CalcMapping(mapping,chrom,start,end,coverage=self.maxHeight,flag=self.flag, threads=self.threads)
        self.CalcPlot=CalcPlot(mapping,chrom,start,end,coverage=self.maxHeight,flag=self.flag, threads=self.threads)
        self.direction=direction
        self.schematic=schematic

    def CalcFRWD(self):
        span=self.end-self.start
        multi=[]
        start=self.start
        end=self.start+self.chunksize
        for _ in range(0,int(span/self.chunksize)):
            multi.append((self.chrom,start,end))
            start=end
            end=end+self.chunksize
        if span/self.chunksize - int(span/self.chunksize)!=0:
            multi.append((self.chrom,start,int(start+((span/self.chunksize - int(span/self.chunksize))*self.chunksize))))
        with Pool(processes=self.threads) as pool:
            results = pool.starmap(self.CalcMapping.plotListFRWRD, multi)
        return results,multi

    def CalcRVRS(self):
        span=self.end-self.start
        multi=[]
        start=self.start
        end=self.start+self.chunksize
        for _ in range(0,int(span/self.chunksize)):
            multi.append((self.chrom,start,end))
            start=end
            end=end+self.chunksize
        if span/self.chunksize - int(span/self.chunksize)!=0:
            multi.append((self.chrom,start,int(start+((span/self.chunksize - int(span/self.chunksize))*self.chunksize))))
        with Pool(processes=self.threads) as pool:
            results = pool.starmap(self.CalcMapping.plotListRVRS, multi)
        return results,multi

    def Calc(self):
        span=self.end-self.start
        multi=[]
        start=self.start
        end=self.start+self.chunksize
        for _ in range(0,int(span/self.chunksize)):
            multi.append((self.chrom,start,end))
            start=end
            end=end+self.chunksize
        if span/self.chunksize - int(span/self.chunksize)!=0:
            multi.append((self.chrom,start,int(start+((span/self.chunksize - int(span/self.chunksize))*self.chunksize))))
        with Pool(processes=self.threads) as pool:
            results = pool.starmap(self.CalcMapping.plotList, multi)
        return results,multi




    def PlotSchematicDir(self):
        resultsR,multiR=self.CalcRVRS()
        resultsF,multiF=self.CalcFRWD()

        self.CalcPlot.startPlot(len(resultsR),self.direction)

        for p,chunk in enumerate(zip(resultsR,multiR)):

            if resultsR[p]==[] and len(resultsR)==1:
                self.CalcPlot.plotEmptyChunk(self.CalcPlot.ax[0,p],chunk[1][1],chunk[1][2])
                self.CalcPlot.ax[0,p].get_xaxis().set_visible(False)
                continue

            if resultsR[p]==[]:
                self.CalcPlot.plotEmptyChunk(self.CalcPlot.ax[0,p],chunk[1][1],chunk[1][2])
                self.CalcPlot.ax[0,p].get_xaxis().set_visible(False)



            if len(resultsR)==1:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                self.CalcPlot.plotChunk(d,self.CalcPlot.ax[0],chunk[1][1],chunk[1][2])
                self.CalcPlot.ax[0].get_xaxis().set_visible(False)
                continue

            if p==0:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                self.CalcPlot.plotChunk(d,self.CalcPlot.ax[0,p],chunk[1][1],chunk[1][2])
                self.CalcPlot.ax[0,p].get_xaxis().set_visible(False)

            else:
                if resultsR[p-1]==[]:
                    d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                    self.CalcPlot.plotChunk(d,self.CalcPlot.ax[0,p],chunk[1][1],chunk[1][2])
                    self.CalcPlot.ax[0,p].get_xaxis().set_visible(False)

                else:
                    d=self.CalcMapping.PlotreadsDF(resultsR[p-1],resultsR[p],multiR[p-1][2])
                    self.CalcPlot.plotChunk(d,self.CalcPlot.ax[0,p],chunk[1][1],chunk[1][2])
                    self.CalcPlot.ax[0,p].get_xaxis().set_visible(False)

        for p,chunk in enumerate(zip(resultsF,multiF)):

            if resultsF[p]==[] and len(resultsR)==1:
                self.CalcPlot.plotEmptyChunk(self.CalcPlot.ax[1,0],chunk[1][1],chunk[1][2])
                continue

            if resultsF[p]==[]:
                self.CalcPlot.plotEmptyChunk(self.CalcPlot.ax[1,p],chunk[1][1],chunk[1][2])


            if len(resultsF)==1:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                self.CalcPlot.plotChunk(d,self.CalcPlot.ax[1],chunk[1][1],chunk[1][2])
                continue

            if p==0:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                self.CalcPlot.plotChunk(d,self.CalcPlot.ax[1,p],chunk[1][1],chunk[1][2])
            else:
                if resultsF[p-1]==[]:
                    d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                    self.CalcPlot.plotChunk(d,self.CalcPlot.ax[1,p],chunk[1][1],chunk[1][2])
                else:
                    d=self.CalcMapping.PlotreadsDF(resultsF[p-1],resultsF[p],multiF[p-1][2])
                    self.CalcPlot.plotChunk(d,self.CalcPlot.ax[1,p],chunk[1][1],chunk[1][2])

        i=self.mapping.split('/')[-1].split('Aligned')[0]
        plt.savefig('{}_{}_{}_{}.pdf'.format(i,self.chrom,str(self.start),str(self.end)))
        plt.close()


    def PlotNucDir(self):
        resultsR,multiR=self.CalcRVRS()
        resultsF,multiF=self.CalcFRWD()
        if len(multiR)>1:
            print('nuc plots over 1000nt not supported due to fontsize')
            sys.exit()
        self.CalcPlot.startPlot(len(resultsR),self.direction)

        for p,chunk in enumerate(zip(resultsR,multiR)):


            if resultsR[p]==[] and len(resultsR)==1:
                self.CalcPlot.plotEmptyChunk(self.CalcPlot.ax[0,p],chunk[1][1],chunk[1][2])
                self.CalcPlot.ax[0,p].get_xaxis().set_visible(False)
                continue

            if resultsR[p]==[]:
                self.CalcPlot.plotEmptyChunk(self.CalcPlot.ax[0,p],chunk[1][1],chunk[1][2])
                self.CalcPlot.ax[0,p].get_xaxis().set_visible(False)
                continue



            if len(resultsR)==1:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                self.CalcPlot.plotChunk(d,self.CalcPlot.ax[0],chunk[1][1],chunk[1][2])
                self.CalcPlot.PlotNucChunk(d,self.CalcPlot.ax[0],chunk[1][1],chunk[1][2],flag=self.flag)
                self.CalcPlot.ax[0].get_xaxis().set_visible(False)
                continue


        for p,chunk in enumerate(zip(resultsF,multiF)):

            if resultsF[p]==[] and len(resultsR)==1:
                self.CalcPlot.plotEmptyChunk(self.CalcPlot.ax[1,0],chunk[1][1],chunk[1][2])
                continue

            if resultsF[p]==[]:
                self.CalcPlot.plotEmptyChunk(self.CalcPlot.ax[1,p],chunk[1][1],chunk[1][2])


            if len(resultsF)==1:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                self.CalcPlot.plotChunk(d,self.CalcPlot.ax[1],chunk[1][1],chunk[1][2])
                self.CalcPlot.PlotNucChunk(d,self.CalcPlot.ax[1],chunk[1][1],chunk[1][2],flag=self.flag)
                continue


        i=self.mapping.split('/')[-1].split('Aligned')[0]
        plt.savefig('{}_{}_{}_{}.pdf'.format(i,self.chrom,str(self.start),str(self.end)))
        plt.close()


    def Plot(self):
        results,multi=self.Calc()
        self.CalcPlot.startPlot(len(results))
        for p,chunk in enumerate(zip(results,multi)):

            if results[p]==[] and len(results)==1:
                self.CalcPlot.plotEmptyChunk(self.CalcPlot.ax,chunk[1][1],chunk[1][2])
                continue

            if results[p]==[]:
                self.CalcPlot.plotEmptyChunk(self.CalcPlot.ax[p],chunk[1][1],chunk[1][2])


            if len(results)==1:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                self.CalcPlot.plotChunk(d,self.CalcPlot.ax,chunk[1][1],chunk[1][2])
                self.CalcPlot.PlotNucChunk(d,self.CalcPlot.ax,chunk[1][1],chunk[1][2],flag=self.flag)
                continue

            if p==0:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                self.CalcPlot.plotChunk(d,self.CalcPlot.ax[p],chunk[1][1],chunk[1][2])
            else:
                if results[p-1]==[]:
                    d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                    self.CalcPlot.plotChunk(d,self.CalcPlot.ax[p],chunk[1][1],chunk[1][2])
                else:
                    d=self.CalcMapping.PlotreadsDF(results[p-1],results[p],multi[p-1][2])
                    self.CalcPlot.plotChunk(d,self.CalcPlot.ax[p],chunk[1][1],chunk[1][2])
        i=self.mapping.split('/')[-1].split('Aligned')[0]
        plt.savefig('{}_{}_{}_{}.pdf'.format(i,self.chrom,str(self.start),str(self.end)))
        plt.close()
