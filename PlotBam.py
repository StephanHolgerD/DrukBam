import argparse
import sys
from classes.MapPlot import PlotMapping

parser = argparse.ArgumentParser(description='Quick reference free plotting of bam file')
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-b','--bam', help='Pos. sorted and indexed bam file', required=True)
requiredNamed.add_argument('-c','--chromosome', help='name of chromosome/contig', required=True)
requiredNamed.add_argument('-s','--start' ,help='start of region of interest', required=True)
requiredNamed.add_argument('-e', '--end',help='end of the region of interest', required=True)


optArguments = parser.add_argument_group('optional arguments')
optArguments.add_argument('--threads',default=1, help="number of cpu's  to run in paralell, ROI <1000 will always use 1 core")
optArguments.add_argument('--direction',default=False,action='store_true', help="split reads by forward and reverse")
optArguments.add_argument('--schematic',default=False,action='store_true', help='plot no nucleotide, recommended for ROI>1000')
optArguments.add_argument('--style',default='classic', help='different style options for the plot, def:classic',choices=['classic','fancy'])
optArguments.add_argument('--flag',default='None', help='different style options for the plot, def:classic',choices=['None','MateUnmapped','SoftClipped','MateUnmappedSoftClipped','Insertion'])




if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args=parser.parse_args()
ploter=PlotMapping(args.bam,args.chromosome,int(args.start),int(args.end),flag=args.flag,schematic=args.schematic,direction=args.direction)
ploter.PlotNucDir()
