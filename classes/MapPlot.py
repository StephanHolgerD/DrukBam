import pysam
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from multiprocessing import Pool
from classes.bamCalc import CalcMapping
from classes.PlotCalc import CalcPlot
import sys
class PlotMapping():
    def __init__(self,mapping,chrom,start,end,coverage=200,flag='None',direction=False,schematic=False,chunksize=1000,threads=1):
        self.mapping=mapping
        self.chrom=chrom
        self.start=start-1
        self.end=end
        self.Fontsize=3
        self.threads=threads
        self.flag=flag
        self.chunksize=chunksize
        self.direction=direction
        if direction:
            self.maxHeight=coverage
        else:
            self.maxHeight=coverage
        self.schematic=schematic
        self.threads=threads
        self.CalcMapping=CalcMapping(mapping,chrom,start,end,coverage=self.maxHeight,flag=self.flag, threads=self.threads)
        self.CalcPlot=CalcPlot(mapping,chrom,start,end,coverage=self.maxHeight,flag=self.flag, threads=self.threads)

    def Plot(self):
        if self.schematic and not self.direction:
            self.PlotSchmematic()
            sys.exit()
        if self.schematic and self.direction:
            self.PlotSchematicDir()
            sys.exit()
        if self.direction and not self.schematic:
            self.PlotNucDir()
            sys.exit()
        if not self.direction and not self.schematic:
            self.PlotNuc()
            sys.exit()



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

        self.CalcPlot.startPlot(len(resultsR),self.direction,self.schematic)

        for p,chunk in enumerate(zip(resultsR,multiR)):
            direction='Reverse'
            if resultsR[p]==[] and len(resultsR)==1:
                self.CalcPlot.AxSet(self.CalcPlot.ax[0,p],chunk[1][1],chunk[1][2],direction=direction)
                self.CalcPlot.ax[0,p].spines['top'].set_visible(False)
                self.CalcPlot.ax[0,p].get_xaxis().set_visible(False)
                continue

            if len(resultsR)==1:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                self.CalcPlot.AxSet(self.CalcPlot.ax[0],chunk[1][1],chunk[1][2],direction=direction)
                self.CalcPlot.plotChunk(d,self.CalcPlot.ax[0],chunk[1][1],chunk[1][2])
                self.CalcPlot.ax[0].spines['top'].set_visible(False)
                self.CalcPlot.ax[0].get_xaxis().set_visible(False)
                continue

            if p==0:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                self.CalcPlot.AxSet(self.CalcPlot.ax[0,p],chunk[1][1],chunk[1][2],direction=direction)
                self.CalcPlot.plotChunk(d,self.CalcPlot.ax[0,p],chunk[1][1],chunk[1][2])
                self.CalcPlot.ax[0,p].get_xaxis().set_visible(False)
                self.CalcPlot.ax[0,p].spines['top'].set_visible(False)
                self.CalcPlot.ax[0,p].get_xaxis().set_visible(False)
                continue

            else:
                if resultsR[p-1]==[]:
                    d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                    self.CalcPlot.AxSet(self.CalcPlot.ax[0,p],chunk[1][1],chunk[1][2],chunk=True)
                    self.CalcPlot.plotChunk(d,self.CalcPlot.ax[0,p],chunk[1][1],chunk[1][2])
                    self.CalcPlot.ax[0,p].spines['top'].set_visible(False)
                    self.CalcPlot.ax[0,p].get_xaxis().set_visible(False)

                else:
                    d=self.CalcMapping.PlotreadsDF(resultsR[p-1],resultsR[p],multiR[p-1][2])
                    self.CalcPlot.AxSet(self.CalcPlot.ax[0,p],chunk[1][1],chunk[1][2],chunk=True)
                    self.CalcPlot.plotChunk(d,self.CalcPlot.ax[0,p],chunk[1][1],chunk[1][2])
                    self.CalcPlot.ax[0,p].spines['top'].set_visible(False)
                    self.CalcPlot.ax[0,p].get_xaxis().set_visible(False)

        for p,chunk in enumerate(zip(resultsF,multiF)):
            direction='Forward'


            if resultsF[p]==[] and len(resultsR)==1:
                self.CalcPlot.AxSet(self.CalcPlot.ax[1,p],chunk[1][1],chunk[1][2],direction=direction)
                continue

            if resultsF[p]==[]:
                self.CalcPlot.AxSet(self.CalcPlot.ax[1,p],chunk[1][1],chunk[1][2],direction=direction)


            if len(resultsF)==1:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                self.CalcPlot.AxSet(self.CalcPlot.ax[1],chunk[1][1],chunk[1][2],direction=direction)
                self.CalcPlot.ax[1].spines['bottom'].set_visible(False)
                self.CalcPlot.ax[1].spines['top'].set_visible(False)
                self.CalcPlot.plotChunk(d,self.CalcPlot.ax[1],chunk[1][1],chunk[1][2])
                continue

            if p==0:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                self.CalcPlot.AxSet(self.CalcPlot.ax[1,p],chunk[1][1],chunk[1][2],direction=direction,chunk=True)
                self.CalcPlot.ax[1,p].spines['bottom'].set_visible(False)
                self.CalcPlot.ax[1,p].spines['top'].set_visible(False)
                self.CalcPlot.plotChunk(d,self.CalcPlot.ax[1,p],chunk[1][1],chunk[1][2])
            else:
                if resultsF[p-1]==[]:
                    d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                    self.CalcPlot.AxSet(self.CalcPlot.ax[1,p],chunk[1][1],chunk[1][2],chunk=True)
                    self.CalcPlot.ax[1,p].spines['bottom'].set_visible(False)
                    self.CalcPlot.ax[1,p].spines['top'].set_visible(False)
                    self.CalcPlot.plotChunk(d,self.CalcPlot.ax[1,p],chunk[1][1],chunk[1][2])
                else:
                    d=self.CalcMapping.PlotreadsDF(resultsF[p-1],resultsF[p],multiF[p-1][2])
                    self.CalcPlot.AxSet(self.CalcPlot.ax[1,p],chunk[1][1],chunk[1][2],chunk=True)
                    self.CalcPlot.ax[1,p].spines['bottom'].set_visible(False)
                    self.CalcPlot.ax[1,p].spines['top'].set_visible(False)
                    self.CalcPlot.plotChunk(d,self.CalcPlot.ax[1,p],chunk[1][1],chunk[1][2])

        i=self.mapping.split('/')[-1].split('Aligned')[0]
        plt.savefig('{}_{}_{}_{}.pdf'.format(i,self.chrom,str(self.start),str(self.end)),bbox_inches='tight')
        plt.close()


    def PlotNucDir(self):
        resultsR,multiR=self.CalcRVRS()
        resultsF,multiF=self.CalcFRWD()
        if len(multiR)>1:
            print('nuc plots over 1000nt not supported due to fontsize')
            sys.exit()
        self.CalcPlot.startPlot(len(resultsR),self.direction,self.schematic)

        for p,chunk in enumerate(zip(resultsR,multiR)):
            direction='Reverse'


            if resultsR[p]==[]:
                self.CalcPlot.AxSet(self.CalcPlot.ax[0],chunk[1][1],chunk[1][2], direction=direction)
                self.CalcPlot.ax[0].get_xaxis().set_visible(False)
                self.CalcPlot.ax[0].spines['bottom'].set_visible(False)

                continue

            if len(resultsR)==1:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                self.CalcPlot.AxSet(self.CalcPlot.ax[0],chunk[1][1],chunk[1][2], direction=direction)
                self.CalcPlot.PlotNucChunk(d,self.CalcPlot.ax[0],chunk[1][1],chunk[1][2],flag=self.flag)
                self.CalcPlot.ax[0].get_xaxis().set_visible(False)
                self.CalcPlot.ax[0].spines['bottom'].set_visible(False)

                continue


        for p,chunk in enumerate(zip(resultsF,multiF)):
            direction='Forward'
            if resultsF[p]==[]:
                self.CalcPlot.AxSet(self.CalcPlot.ax[1],chunk[1][1],chunk[1][2], direction=direction)
                continue

            if len(resultsF)==1:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                self.CalcPlot.ax[1].spines['bottom'].set_visible(False)
                self.CalcPlot.AxSet(self.CalcPlot.ax[1],chunk[1][1],chunk[1][2], direction=direction)
                self.CalcPlot.PlotNucChunk(d,self.CalcPlot.ax[1],chunk[1][1],chunk[1][2],flag=self.flag)
                continue


        i=self.mapping.split('/')[-1].split('Aligned')[0]
        plt.savefig('{}_{}_{}_{}.pdf'.format(i,self.chrom,str(self.start),str(self.end)),bbox_inches='tight')
        plt.close()

    def PlotNuc(self):
        results,multi=self.Calc()
        self.CalcPlot.startPlot(len(results),self.direction,self.schematic)
        if len(multi)>1:
            print('nuc plots over 1000nt not supported due to fontsize')
            sys.exit()
        for p,chunk in enumerate(zip(results,multi)):

            if results[p]==[]:

                self.CalcPlot.AxSet(self.CalcPlot.ax,chunk[1][1],chunk[1][2])
                continue
            else:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                self.CalcPlot.AxSet(self.CalcPlot.ax,chunk[1][1],chunk[1][2])
                self.CalcPlot.ax.spines['bottom'].set_visible(False)
                self.CalcPlot.PlotNucChunk(d,self.CalcPlot.ax,chunk[1][1],chunk[1][2],flag=self.flag)
                continue



        i=self.mapping.split('/')[-1].split('Aligned')[0]
        plt.savefig('{}_{}_{}_{}.pdf'.format(i,self.chrom,str(self.start),str(self.end)),bbox_inches='tight')
        plt.close()

    def PlotSchmematic(self):
        results,multi=self.Calc()
        self.CalcPlot.startPlot(len(results),self.direction,self.schematic)
        for p,chunk in enumerate(zip(results,multi)):

            if results[p]==[] and len(results)==1:
                self.CalcPlot.AxSet(self.CalcPlot.ax,chunk[1][1],chunk[1][2])
                continue

            if results[p]==[]:
                self.CalcPlot.AxSet(self.CalcPlot.ax[p],chunk[1][1],chunk[1][2])


            if len(results)==1:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                self.CalcPlot.AxSet(self.CalcPlot.ax,chunk[1][1],chunk[1][2])
                self.CalcPlot.ax.spines['bottom'].set_visible(False)
                self.CalcPlot.plotChunk(d,self.CalcPlot.ax,chunk[1][1],chunk[1][2])

                #self.CalcPlot.PlotNucChunk(d,self.CalcPlot.ax,chunk[1][1],chunk[1][2],flag=self.flag)
                continue

            if p==0:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                self.CalcPlot.AxSet(self.CalcPlot.ax[p],chunk[1][1],chunk[1][2],chunk=True)
                self.CalcPlot.ax[p].spines['bottom'].set_visible(False)
                self.CalcPlot.plotChunk(d,self.CalcPlot.ax[p],chunk[1][1],chunk[1][2])
            else:
                if results[p-1]==[]:
                    d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                    self.CalcPlot.ax[p].spines['bottom'].set_visible(False)
                    self.CalcPlot.AxSet(self.CalcPlot.ax[p],chunk[1][1],chunk[1][2],chunk=True)
                    self.CalcPlot.plotChunk(d,self.CalcPlot.ax[p],chunk[1][1],chunk[1][2])
                else:
                    d=self.CalcMapping.PlotreadsDF(results[p-1],results[p],multi[p-1][2])
                    self.CalcPlot.ax[p].spines['bottom'].set_visible(False)
                    self.CalcPlot.AxSet(self.CalcPlot.ax[p],chunk[1][1],chunk[1][2],chunk=True)
                    self.CalcPlot.plotChunk(d,self.CalcPlot.ax[p],chunk[1][1],chunk[1][2])



        i=self.mapping.split('/')[-1].split('Aligned')[0]
        plt.savefig('{}_{}_{}_{}.pdf'.format(i,self.chrom,str(self.start),str(self.end)),bbox_inches='tight')
        plt.close()
