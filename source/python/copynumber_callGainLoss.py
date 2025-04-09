#!/usr/bin/env python3

"""
Script to use stringent cutoffs to make GAIN or LOSS copy number region calls.
Based on recommendations from Chen et al. (bioRxiv 2021), the default 
gain cutoff is 3 and the loss cutoff is 1.5.

Takes segment file from Sequenza (segments.txt), FACETS (.cncf), or CNVkit (.cns)

Sequenza header: use ("chromosome","start.pos","end.pos","CNt")
"chromosome"    "start.pos"  "end.pos"   "Bf"    "N.BAF" "sd.BAF"    "depth.ratio"   "N.ratio"   "sd.ratio"  "CNt"   "A" "B" "LPP"


FACETS header: use (chrom,start,end,tcn.em)
NOTE: chroms are '1' rather than 'chr1' --convert to chr1
chrom   seg num.mark    nhet    cnlr.median mafR    segclust    cnlr.median.clust   mafR.clust  start   end cf.em   tcn.em  lcn.em


CNVkit header: use (chromosome,start,end,cn)
chromosome  start   end gene    log2    baf cn  cn1 cn2 depth   p_ttest probes  weight
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f [segmentation file from either Sequenza, Facets, or CNVkit]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--file", help="segment file from either Sequenza, Facets, or CNVkit file")
    optparser.add_option("-o", "--output", help="output gain/loss calls with stringent cutoffs")
    optparser.add_option("-g", "--gain_cutoff", help="total copy number cutoff for gains, default: 3", default="3")
    optparser.add_option("-l", "--loss_cutoff", help="total copy number cutoff for losses, default: 1.5", default="1.5")

    (options, args) = optparser.parse_args(sys.argv)

    if not options.file or not options.output:
        optparser.print_help()
        sys.exit(-1)

    gain_cut = float(options.gain_cutoff)
    loss_cut = float(options.loss_cutoff)

    with open(options.file) as f, open(options.output,"w") as out:
        hdr = f.readline().strip().split("\t")
        ## Write header for output
        out.write("%s\n" % "\t".join(["chrom","start","end","total_cn","call"]))
        if 'tcn.em' in hdr: ## FACETS file
            caller = "facets"
            chr_i = 'chrom'
            st_i = 'start'
            end_i = 'end'
            tcn_i = 'tcn.em'
        elif '\"CNt\"' in hdr: ## Sequenza file
            caller = "sequenza"
            chr_i = '\"chromosome\"'
            st_i = '\"start.pos\"'
            end_i = '\"end.pos\"'
            tcn_i = '\"CNt\"'
        else: ## CNVKit file
            caller = "cnvkit"
            chr_i = 'chromosome'
            st_i = 'start'
            end_i = 'end'
            tcn_i = 'cn'


        ## Read the file
        for l in f:
            tmp = dict(zip(hdr, l.strip().split("\t")))
            if tmp[tcn_i] in ["NA", ""]:  
              continue  
            copyNum = int(float(tmp[tcn_i]))
            ## Check cutoffs
            if copyNum >= gain_cut or copyNum <= loss_cut:
                ## Handle chrom name oddities
                if caller == "facets":
                    chrom = "chr%s" % tmp[chr_i] if tmp[chr_i] != '23' else 'chrX'
                elif caller == "sequenza": ## Sequenza remove quoutes eg. "chr1"
                    chrom = eval(tmp[chr_i])
                else:
                    chrom = tmp[chr_i]
                call = "GAIN" if copyNum >= gain_cut else "LOSS"

                start = tmp[st_i]
                end = tmp[end_i]

                ## Write to output
                out.write("%s\n" % "\t".join([chrom,start,end,str(copyNum),call]))

if __name__=='__main__':
    main()
