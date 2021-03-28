import pysam
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from multiprocessing import Pool
from DrukBam.bamCalc import CalcMapping
from DrukBam.PlotCalc import CalcPlot
import sys
class PlotMapping():
    def __init__(self,mapping,chrom,start,end,coverage=200,flag='None',direction=False,schematic=False,chunksize=1000,threads=1,fasta=None,output='current working directory',
    out_name='name of mapping',vcf=False,style='classic',outfmt='pdf',outlineoff=False):
        self.outputdir=output
        self.out_name=out_name
        self.mapping=mapping
        self.chrom=chrom
        self.start=start
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
        self.fasta=fasta
        self.style=style
        self.CalcMapping=CalcMapping(mapping,chrom,start,end,coverage=self.maxHeight,flag=self.flag, threads=self.threads)
        self.CalcPlot=CalcPlot(mapping,chrom,self.start,self.end,coverage=self.maxHeight,flag=self.flag, threads=self.threads,fasta=self.fasta,style=self.style,outlineoff=outlineoff)
        self.vcf=vcf
        self.outfmt=outfmt
    def savePlot(self,schmematic=False,direction=False,reference=False):
        if self.out_name=='name of mapping':
            o=self.mapping.split('/')[-1]
        else:
            o=self.out_name
        if self.outputdir=='current working directory':
            od='./'
        else:
            od=self.outputdir
        r=''
        d=''
        s=''
        if schmematic:
            r='S'
        if direction:
            d='D'
        if reference != 'None':
            r='R'

        self.fig.savefig('{}/{}_{}_{}_{}{}{}{}.{}'.format(od,o,self.chrom,str(self.start),str(self.end),s,d,r,self.outfmt),bbox_inches='tight',dpi=400)
        plt.close()

    def Plot(self):
        if self.schematic and not self.direction:
            self.PlotSchmematic()
            self.savePlot(schmematic=True)
            return
        if self.schematic and self.direction:
            self.PlotSchematicDir()
            self.savePlot(schmematic=True,direction=True)
            return
        if self.direction and not self.schematic:
            if (self.end - self.start) > self.chunksize:
                print('span larger than chunksize, please increase chunksize, calc. speed will slow down with large chunksizes')
                sys.exit()
            self.PlotNucDir()
            self.savePlot(direction=True,reference=self.fasta)
            return
        if not self.direction and not self.schematic:
            if (self.end - self.start) > self.chunksize:
                print('span larger than chunksize, please increase chunksize, calc. speed will slow down with large chunksizes')
                sys.exit()
            self.PlotNuc()
            self.savePlot(reference=self.fasta)
            return


    def CalcFRWD(self):
        span=self.end-self.start
        multi=[]
        start=self.start
        end=self.start+self.chunksize
        if span<self.chunksize:
            results=[self.CalcMapping.plotListFRWRD(self.chrom,self.start,self.end)]
            multi.append((self.chrom,self.start,self.end))
        else:
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
        if span<self.chunksize:
            results=[self.CalcMapping.plotListRVRS(self.chrom,self.start,self.end)]
            multi.append((self.chrom,self.start,self.end))
        else:
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
        if span < self.chunksize:
            results=[self.CalcMapping.plotList(self.chrom,self.start,self.end)]
            multi.append((self.chrom,self.start,self.end))
        else:
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

        self.fig,self.ax=self.CalcPlot.startPlot(len(resultsR),self.direction,self.schematic)

        for p,chunk in enumerate(zip(resultsR,multiR)):
            direction='Reverse'
            if resultsR[p]==[] and len(resultsR)==1:
                self.CalcPlot.AxSet(self.ax[0,p],chunk[1][1],chunk[1][2],direction=direction)
                self.ax[0,p].spines['top'].set_visible(False)
                self.ax[0,p].get_xaxis().set_visible(False)
                continue

            if len(resultsR)==1:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                self.CalcPlot.AxSet(self.ax[0],chunk[1][1],chunk[1][2],direction=direction)
                self.CalcPlot.plotChunk(d,self.ax[0],chunk[1][1],chunk[1][2])
                self.ax[0].spines['top'].set_visible(False)
                self.ax[0].get_xaxis().set_visible(False)
                continue

            if p==0:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                self.CalcPlot.AxSet(self.ax[0,p],chunk[1][1],chunk[1][2],direction=direction)
                self.CalcPlot.plotChunk(d,self.ax[0,p],chunk[1][1],chunk[1][2])
                self.ax[0,p].get_xaxis().set_visible(False)
                self.ax[0,p].spines['top'].set_visible(False)
                self.ax[0,p].get_xaxis().set_visible(False)
                continue

            else:
                if resultsR[p-1]==[]:
                    d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                    self.CalcPlot.AxSet(self.ax[0,p],chunk[1][1],chunk[1][2],chunk=True)
                    self.CalcPlot.plotChunk(d,self.ax[0,p],chunk[1][1],chunk[1][2])
                    self.ax[0,p].spines['top'].set_visible(False)
                    self.ax[0,p].get_xaxis().set_visible(False)

                else:
                    d=self.CalcMapping.PlotreadsDF(resultsR[p-1],resultsR[p],multiR[p-1][2])
                    self.CalcPlot.AxSet(self.ax[0,p],chunk[1][1],chunk[1][2],chunk=True)
                    self.CalcPlot.plotChunk(d,self.ax[0,p],chunk[1][1],chunk[1][2])
                    self.ax[0,p].spines['top'].set_visible(False)
                    self.ax[0,p].get_xaxis().set_visible(False)

        for p,chunk in enumerate(zip(resultsF,multiF)):
            direction='Forward'


            if resultsF[p]==[] and len(resultsR)==1:
                self.CalcPlot.AxSet(self.ax[1,p],chunk[1][1],chunk[1][2],direction=direction)
                continue

            if resultsF[p]==[]:
                self.CalcPlot.AxSet(self.ax[1,p],chunk[1][1],chunk[1][2],direction=direction)


            if len(resultsF)==1:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                self.CalcPlot.AxSet(self.ax[1],chunk[1][1],chunk[1][2],direction=direction)
                self.ax[1].spines['bottom'].set_visible(False)
                self.ax[1].spines['top'].set_visible(False)
                self.CalcPlot.plotChunk(d,self.ax[1],chunk[1][1],chunk[1][2])
                continue

            if p==0:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                self.CalcPlot.AxSet(self.ax[1,p],chunk[1][1],chunk[1][2],direction=direction,chunk=True)
                self.ax[1,p].spines['bottom'].set_visible(False)
                self.ax[1,p].spines['top'].set_visible(False)
                self.CalcPlot.plotChunk(d,self.ax[1,p],chunk[1][1],chunk[1][2])
            else:
                if resultsF[p-1]==[]:
                    d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                    self.CalcPlot.AxSet(self.ax[1,p],chunk[1][1],chunk[1][2],chunk=True)
                    self.ax[1,p].spines['bottom'].set_visible(False)
                    self.ax[1,p].spines['top'].set_visible(False)
                    self.CalcPlot.plotChunk(d,self.ax[1,p],chunk[1][1],chunk[1][2])
                else:
                    d=self.CalcMapping.PlotreadsDF(resultsF[p-1],resultsF[p],multiF[p-1][2])
                    self.CalcPlot.AxSet(self.ax[1,p],chunk[1][1],chunk[1][2],chunk=True)
                    self.ax[1,p].spines['bottom'].set_visible(False)
                    self.ax[1,p].spines['top'].set_visible(False)
                    self.CalcPlot.plotChunk(d,self.ax[1,p],chunk[1][1],chunk[1][2])




    def PlotNucDir(self):
        resultsR,multiR=self.CalcRVRS()
        resultsF,multiF=self.CalcFRWD()
        if len(multiR)>1:
            print('nuc plots over 1000nt not supported due to fontsize')
            sys.exit()
        self.fig,self.ax=self.CalcPlot.startPlot(len(resultsR),self.direction,self.schematic)

        for p,chunk in enumerate(zip(resultsR,multiR)):
            direction='Reverse'


            if resultsR[p]==[]:
                self.CalcPlot.AxSet(self.ax[0],chunk[1][1],chunk[1][2], direction=direction,vcf=self.vcf)
                self.ax[0].get_xaxis().set_visible(False)
                self.ax[0].spines['bottom'].set_visible(False)

                continue

            if len(resultsR)==1:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                self.CalcPlot.AxSet(self.ax[0],chunk[1][1],chunk[1][2], direction=direction,vcf=self.vcf)

                self.CalcPlot.PlotNucChunk(d,self.ax[0],chunk[1][1],chunk[1][2],flag=self.flag)
                self.ax[0].get_xaxis().set_visible(False)
                self.ax[0].spines['bottom'].set_visible(False)

                continue


        for p,chunk in enumerate(zip(resultsF,multiF)):
            direction='Forward'
            if resultsF[p]==[]:
                self.CalcPlot.AxSet(self.ax[1],chunk[1][1],chunk[1][2], direction=direction,vcf=self.vcf)
                continue

            if len(resultsF)==1:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                self.ax[1].spines['bottom'].set_visible(False)
                self.CalcPlot.AxSet(self.ax[1],chunk[1][1],chunk[1][2], direction=direction,vcf=self.vcf)

                self.CalcPlot.PlotNucChunk(d,self.ax[1],chunk[1][1],chunk[1][2],flag=self.flag)

                continue

        self.CalcPlot.PlotFasta(self.ax[0])
        self.CalcPlot.PlotFasta(self.ax[1])


    def PlotNuc(self):
        results,multi=self.Calc()
        self.fig,self.ax=self.CalcPlot.startPlot(len(results),self.direction,self.schematic)
        if len(multi)>1:
            print('nuc plots over 1000nt not supported due to fontsize')
            sys.exit()
        for p,chunk in enumerate(zip(results,multi)):

            if results[p]==[]:

                self.CalcPlot.AxSet(self.ax,chunk[1][1],chunk[1][2])
                continue
            else:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                self.CalcPlot.AxSet(self.ax,chunk[1][1],chunk[1][2],vcf=self.vcf)
                self.ax.spines['bottom'].set_visible(False)
                self.CalcPlot.PlotNucChunk(d,self.ax,chunk[1][1],chunk[1][2],flag=self.flag)
                continue
        self.CalcPlot.PlotFasta(self.ax)





    def PlotSchmematic(self):
        results,multi=self.Calc()
        self.fig,self.ax=self.CalcPlot.startPlot(len(results),self.direction,self.schematic)
        for p,chunk in enumerate(zip(results,multi)):

            if results[p]==[] and len(results)==1:
                self.CalcPlot.AxSet(self.ax,chunk[1][1],chunk[1][2])
                continue

            if results[p]==[]:
                self.CalcPlot.AxSet(self.ax[p],chunk[1][1],chunk[1][2])


            if len(results)==1:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                self.CalcPlot.AxSet(self.ax,chunk[1][1],chunk[1][2])
                self.ax.spines['bottom'].set_visible(False)
                self.CalcPlot.plotChunk(d,self.ax,chunk[1][1],chunk[1][2])

                #self.CalcPlot.PlotNucChunk(d,self.ax,chunk[1][1],chunk[1][2],flag=self.flag)
                continue

            if p==0:
                d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                self.CalcPlot.AxSet(self.ax[p],chunk[1][1],chunk[1][2],chunk=True)
                self.ax[p].spines['bottom'].set_visible(False)
                self.CalcPlot.plotChunk(d,self.ax[p],chunk[1][1],chunk[1][2])
            else:
                if results[p-1]==[]:
                    d=pd.DataFrame(chunk[0],columns=['y','start','end','direction','name','qSeq','cigar','mateMap'])
                    self.ax[p].spines['bottom'].set_visible(False)
                    self.CalcPlot.AxSet(self.ax[p],chunk[1][1],chunk[1][2],chunk=True)
                    self.CalcPlot.plotChunk(d,self.ax[p],chunk[1][1],chunk[1][2])
                else:
                    d=self.CalcMapping.PlotreadsDF(results[p-1],results[p],multi[p-1][2])
                    self.ax[p].spines['bottom'].set_visible(False)
                    self.CalcPlot.AxSet(self.ax[p],chunk[1][1],chunk[1][2],chunk=True)
                    self.CalcPlot.plotChunk(d,self.ax[p],chunk[1][1],chunk[1][2])
