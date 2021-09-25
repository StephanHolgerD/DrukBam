#!/usr/bin/env python3
import argparse
import sys
from DrukBam.MapPlot import PlotMapping
from DrukBam.vcfParse import VcfPlotter
from time import time
#from importlib import resources
def main():
    parser = argparse.ArgumentParser(description='Quick reference free plotting of bam file')
    subparsers = parser.add_subparsers(help='plot variants from vcf file with padding / plot region')

    vcf_parser = subparsers.add_parser("vcf")
    region_parser = subparsers.add_parser("region")

    requiredNamed = region_parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-b','--bam', help='Pos. sorted and indexed bam file', required=True)
    requiredNamed.add_argument('-c','--chromosome', help='name of chromosome/contig', required=True)
    requiredNamed.add_argument('-s','--start' ,help='start of region of interest', required=True)
    requiredNamed.add_argument('-e', '--end',help='end of the region of interest', required=True)

    requiredNamed = vcf_parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-b','--bam', help='Pos. sorted and indexed bam file', required=True)
    requiredNamed.add_argument('-v','--vcf', help='vcf file with variants of interest', required=True)
    requiredNamed.add_argument('-p','--padding', help='number of nt around the variant', default=25,type=int)
    requiredNamed.add_argument('--highlight', help='highlight the position of interest', default=True, action='store_true')


    optArguments = region_parser.add_argument_group('optional arguments')
    optArguments.add_argument('--threads',default=1, help="number of cpu's  to run in paralell, ROI <1000 will always use 1 core",type=int)
    optArguments.add_argument('--maxcoverage',default=200, help="max cov to plot",type=int)

    optArguments.add_argument('--direction',default=False,action='store_true', help="split reads by forward and reverse")
    optArguments.add_argument('--schematic',default=False,action='store_true', help='plot no nucleotide, recommended for ROI>1000')
    optArguments.add_argument('--style',default=None, help='different style options for the plot, provide .ini file')
    optArguments.add_argument('--Flag',default=None, help='provide flag of non-interest, all single flags are not allowed to be in total flag of read')
    optArguments.add_argument('--flag',default=None, help='provide flag of interest, all single flags needs to be in toctal flag of read')

    optArguments.add_argument('--fasta',default='None', help='fasta file for reference related plotting')
    optArguments.add_argument('--outputdir',default='current working directory', help='directory for output')
    optArguments.add_argument('-i','--id',default='name of mapping', help='output filename')
    optArguments.add_argument('--chunksize',default=1000, help='max size of visualized area, can be increases but will sow down calculation',type=int)
    optArguments.add_argument('--outfmt',default='pdf', help='format of plot, choose between pdf,svg,png')
    optArguments.add_argument('--outlineoff',default=False, help='plotting of read outline',action='store_true')




    optArguments = vcf_parser.add_argument_group('optional arguments')
    optArguments.add_argument('--threads',default=1, help="number of cpu's  to run in paralell, ROI <1000 will always use 1 core",type=int)
    optArguments.add_argument('--maxcoverage',default=200, help="max cov to plot",type=int)

    optArguments.add_argument('--direction',default=False,action='store_true', help="split reads by forward and reverse")
    optArguments.add_argument('--schematic',default=False,action='store_true', help='plot no nucleotide, recommended for ROI>1000')
    optArguments.add_argument('--style',default=None, help='different style options for the plot, provide .ini file')
    optArguments.add_argument('--flag',default=None, help='provide flag of interest, all single flags needs to be in toctal flag of read')
    optArguments.add_argument('--Flag',default=None, help='provide flag of non-interest, all single flags are not allowed to be in total flag of read')

    optArguments.add_argument('--fasta',default='None', help='fasta file for reference related plotting')
    optArguments.add_argument('--outputdir',default='current working directory', help='directory for output')
    optArguments.add_argument('-i','--id',default='name of mapping', help='output filename')
    optArguments.add_argument('--chunksize',default=1000, help='max size of visualized area, can be increases but will sow down calculation',type=int)
    optArguments.add_argument('--outfmt',default='pdf', help='format of plot, choose between pdf,svg,png')
    optArguments.add_argument('--outlineoff',default=False, help='plotting of read outline',action='store_true')


    t=time()
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args=parser.parse_args()
    if args.__contains__('start'):
        ploter=PlotMapping(args.bam,args.chromosome,
        int(args.start),
        int(args.end),
        schematic=args.schematic,
        direction=args.direction,
        coverage=args.maxcoverage,threads=args.threads,
        fasta=args.fasta,
        out_name=args.id,
        output=args.outputdir,
        chunksize=args.chunksize,
        style=args.style,
        outfmt=args.outfmt,
        outlineoff=args.outlineoff,
        flag=args.flag,
        Flag=args.Flag)
        ploter.Plot()
    if args.__contains__('vcf'):
        ploter=VcfPlotter(args.vcf,args.bam,
        schematic=args.schematic,
        direction=args.direction,
        coverage=args.maxcoverage,threads=args.threads,
        fasta=args.fasta,
        out_name=args.id,
        output=args.outputdir,
        padding=args.padding,
        chunksize=args.chunksize,
        style=args.style,
        outfmt=args.outfmt,
        outlineoff=args.outlineoff,
        Flag=args.Flag,
        flag=args.flag)
        ploter.MultiPlot()
        print(time()-t)

if __name__ == "__main__":
    main()
