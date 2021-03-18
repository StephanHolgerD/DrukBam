import pysam
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from multiprocessing import Pool
from druk_bam.MapPlot import PlotMapping
import sys
from time import time
class VcfPlotter():
    def __init__(self,vcf,mapping,coverage=200,flag='None',padding=20,direction=False,schematic=False,threads=1,fasta=None,output='current working directory',out_name='name of mapping'):
        self.outputdir=output
        self.out_name=out_name
        self.mapping=mapping
        self.vcf=vcf
        self.padding=padding
        self.Fontsize=3
        self.threads=threads
        self.flag=flag
        self.direction=direction
        if direction:
            self.maxHeight=coverage
        else:
            self.maxHeight=coverage
        self.schematic=schematic
        self.threads=threads
        self.fasta=fasta
    def Plot(self):
        t=time()
        with pysam.VariantFile(self.vcf) as v:
            for record in v:
                ploter=PlotMapping(
                self.mapping,
                record.chrom,
                int(record.pos-self.padding),
                int(record.pos+self.padding),
                flag=self.flag,
                schematic=self.schematic,
                direction=self.direction,
                coverage=self.maxHeight,
                threads=self.threads,
                fasta=self.fasta,
                out_name=self.out_name,
                output=self.outputdir)
                ploter.Plot()
                print(time()-t)
    def PlotV(self,c,s,e):
            ploter=PlotMapping(
                self.mapping,
                c,
                int(s),
                int(e),
                flag=self.flag,
                schematic=self.schematic,
                direction=self.direction,
                coverage=self.maxHeight,
                threads=self.threads,
                fasta=self.fasta,
                out_name=self.out_name,
                output=self.outputdir)
            ploter.Plot()
    def MultiPlot(self):
        multi=[]
        with pysam.VariantFile(self.vcf) as v:
            for record in v:
                multi.append((record.chrom,(record.pos-self.padding),(record.pos+self.padding)))
        with Pool(processes=self.threads) as pool:
            results = pool.starmap(self.PlotV, multi)
