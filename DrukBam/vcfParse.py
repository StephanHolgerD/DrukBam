import pysam
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from multiprocessing import Pool
from DrukBam.MapPlot import PlotMapping
import sys
from time import time
class VcfPlotter():
    def __init__(self,vcf,mapping,coverage=200,flag=None,Flag=None,chunksize=1000,padding=20,direction=False,schematic=False,threads=1,
    fasta=None,output='current working directory',out_name='name of mapping',style='classic',outfmt='pdf',outlineoff=False,clipped=1):
        self.outputdir=output
        self.out_name=out_name
        self.mapping=mapping
        self.vcf=vcf
        self.padding=padding
        self.Fontsize=3
        self.threads=threads
        self.flag=flag
        self.Flag=Flag
        self.direction=direction
        self.clipped=clipped
        if direction:
            self.maxHeight=coverage
        else:
            self.maxHeight=coverage
        self.schematic=schematic
        self.threads=threads
        self.fasta=fasta
        self.chunksize=chunksize
        self.style=style
        self.outfmt=outfmt
        self.outlineoff=outlineoff
    def PlotV(self,c,s,e):
            ploter=PlotMapping(
                self.mapping,
                c,
                int(s),
                int(e),
                flag=self.flag,
                Flag=self.Flag,
                clipped=self.clipped,
                schematic=self.schematic,
                direction=self.direction,
                coverage=self.maxHeight,
                threads=self.threads,
                fasta=self.fasta,
                out_name=self.out_name,
                output=self.outputdir,
                chunksize=self.chunksize,
                vcf=True,
                outfmt=self.outfmt,
                style=self.style,
                outlineoff=self.outlineoff)
            ploter.Plot()
    def MultiPlot(self):
        multi=[]
        with pysam.VariantFile(self.vcf) as v:
            for record in v:
                multi.append((record.chrom,(record.pos-self.padding),(record.pos+self.padding)))
        print(multi)
        with Pool(processes=self.threads) as pool:
            results = pool.starmap(self.PlotV, multi)
