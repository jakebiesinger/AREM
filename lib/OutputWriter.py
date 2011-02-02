# Time-stamp: <2011-01-20 18:21:42 Jake Biesinger>

"""Description:

Copyright (c) 2008,2009,2010 Yong Zhang, Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).

@status: beta
@version: $Revision$
@originalauthor:  Yong Zhang, Tao Liu
@originalcontact: taoliu@jimmy.harvard.edu

Modifications to probabilistically align reads to regions with highest
enrichment performed by Jacob Biesinger. Repackaged as "AREM" in accordance
with copyright restrictions.

@author: Jacob Biesinger, Daniel Newkirk, Alvin Chon 
@contact: jake.biesinger@gmail.com, dnewkirk@uci.edu, achon@uci.edu

Changes to this file since original release of MACS 1.4 (summer wishes):
  December/January 2011
    * Updated names (AREM, not MACS14)
"""

# ------------------------------------
# python modules
# ------------------------------------
import os
import sys
from array import array

from Constants import *

# ------------------------------------
# Misc functions
# ------------------------------------
def zwig_write (trackI, subdir, fileprefix, d, log=None,space=10, single=False):
    """Write shifted tags information in wiggle file in a given
    step. Then compress it using 'gzip' program.

    trackI: shifted tags from PeakDetect object
    subdir: directory where to put the wiggle file
    fileprefix: wiggle file prefix
    d     : d length
    log   : logging function, default is sys.stderr.write
    space : space to write tag number on spots, default 10
    """
    if not log:
        log = lambda x: sys.stderr.write(x+"\n")
    chrs = trackI.get_chr_names()
    os.makedirs (subdir)
    step = 10000000 + 2*d

    if single:
        log("write to a wiggle file")
        f = os.path.join(subdir,fileprefix+"_all"+".wig")
        wigfhd = open(f,"w")
        wigfhd.write("track type=wiggle_0 name=\"%s_all\" description=\"Shifted Merged AREM tag counts for every %d bp\"\n" % (fileprefix.replace('_afterfiting',''), space)) # data type line        
    
    for chrom in chrs:
        if not single:
            f = os.path.join(subdir,fileprefix+"_"+chrom+".wig")
            log("write to "+f+" for chromosome "+chrom)
            wigfhd = open(f,"w")
            # suggested by dawe
            #wigfhd.write("track type=wiggle_0 name=\"AREM_counts_after_shifting\" description=\"Shifted Merged AREM tag counts for every %d bp\"\n" % (space)) # data type line
            wigfhd.write("track type=wiggle_0 name=\"%s_%s\" description=\"Shifted Merged AREM tag counts for every %d bp\"\n" % (fileprefix.replace('_afterfiting',''), chrom, space)) # data type line
        else:
            log("write data for chromosome "+chrom)
            
        wigfhd.write("variableStep chrom=%s span=%d\n" % (chrom,space))
        tags = trackI.get_locations_by_chr(chrom)[0]
        l = len(tags)
        window_counts = array(BYTE4,[0]*step)
        startp = -1*d
        endp   = startp+step
        index_tag = 0

        while index_tag<l:
            s = tags[index_tag]-d/2     # start of tag
            e = s+d                     # end of tag
            
            if e < endp:
                # project tag to window_counts line
                ps = s-startp # projection start
                pe = ps+d     # projection end
                for i in xrange(ps,pe):
                    window_counts[i] += 1
                index_tag += 1
            else:
                # write it to zwig file then reset parameters
                # keep this tag for next window
                for i in xrange(d,step-d,space):
                    if window_counts[i] == 0:
                        pass
                    else:
                        wigfhd.write("%d\t%d\n" % (i+startp+1,window_counts[i]))
                # reset
                window_counts_next = array(BYTE4,[0]*step)
                # copy d values from the tail of previous window to next window
                for n,i in enumerate(xrange(step-2*d,step)): # debug
                    window_counts_next[n] = window_counts[i]
                window_counts = window_counts_next
                startp = endp - 2*d
                endp = startp+step
        # last window
        for i in xrange(d,step-d,space):
            if window_counts[i] == 0:
                pass
            else:
                wigfhd.write("%d\t%d\n" % (i+startp+1,window_counts[i]))
        if not single:
            wigfhd.close()
            log("compress the wiggle file using gzip...")
            os.system("gzip "+f)
    if single:
        wigfhd.close()
        log("compress the wiggle file using gzip...")
        os.system("gzip "+f)

# def zwig_score_write (treatwig, controlbk, d, log=None, space=10, save_pvalue=True, lambdaset=None, treatnum=10000, controlnum=10000, bglambda=0.01):
#     """Write scores. Then compress it using 'gzip' program.

#     treat: shifted tags from PeakDetect object for treatment
#     control: shifted tags from PeakDetect object for control    
#     subdir: directory where to put the wiggle file
#     fileprefix: wiggle file prefix
#     d     : d length
#     log   : logging function, default is sys.stderr.write
#     space : space to write tag number on spots, default 10
#     pvalue: True: use -10*log(10,pvalue) as score, otherwise use fold-enrichment
#     """
#     if not log:
#         log = lambda x: sys.stderr.write(x+"\n")
#     chrs = treatwig.get_chr_names()
#     (sregion,mregion,lregion)=lambdaset
#     sshift = sregion/2
#     mshift = mregion/2
#     lshift = lregion/2
#     dshift = d/2
#     scores = WigTrackI()
#     scores.span = space
#     space2d = float(space)/d
#     ratio = float(treatnum)/controlnum
#     for chrom in chrs:
    
#         smallregionsum    = 0
#         middleregionsum   = 0
#         largeregionsum    = 0
        
#         (treatpos,treatcount) = treatwig.get_data_by_chr(chrom)
#         controlbk_chr = controlbk[chrom]
#         for i in range(len(treatpos)):
#             pos = treatpos[i]
#             count = treatcount[i]
#             dregion_left = pos-dshift
#             dregion_right = pos+dshift
#             dcontrollambda = sum(controlbk_chr.pp2v(dregion_left,dregion_right))/d*space
#             sregion_left = pos-sshift
#             sregion_right = pos+sshift
#             scontrollambda = sum(controlbk_chr.pp2v(sregion_left,sregion_right))/sregion*space
#             mregion_left = pos-sshift
#             mregion_right = pos+sshift
#             mcontrollambda = sum(controlbk_chr.pp2v(mregion_left,mregion_right))/mregion*space
#             lregion_left = pos-sshift
#             lregion_right = pos+sshift
#             lcontrollambda = sum(controlbk_chr.pp2v(lregion_left,lregion_right))/lregion*space
#             #print pos,count,dcontrollambda,scontrollambda,mcontrollambda,lcontrollambda,bglambda
#             maxlambda = max(dcontrollambda,scontrollambda,mcontrollambda,lcontrollambda,bglambda)
#             p_tmp = poisson_cdf(count,maxlambda*ratio,lower=False)
#             if p_tmp <= 0:
#                 peak_pvalue = 3100
#             else:
#                 peak_pvalue = mathlog(p_tmp,10) * -10
#             if save_pvalue:
#                 scores.add_loc(chrom,pos,peak_pvalue)
#             else:
#                 scores.add_loc(chrom,pos,count/maxlambda)
            
#     return scores

def model2r_script(model,filename,name):
    rfhd = open(filename,"w")
    p = model.plus_line
    m = model.minus_line
    s = model.shifted_line
    d = model.d
    w = len(p)
    norm_p = [0]*w
    norm_m = [0]*w
    norm_s = [0]*w
    sum_p = sum(p)
    sum_m = sum(m)
    sum_s = sum(s)
    for i in range(w):
        norm_p[i] = float(p[i])*100/sum_p
        norm_m[i] = float(m[i])*100/sum_m
        norm_s[i] = float(s[i])*100/sum_s
    rfhd.write("# R script for Peak Model\n")
    rfhd.write("#  -- generated by AREM\n")

    rfhd.write("""p <- c(%s)
m <- c(%s)
s <- c(%s)
x <- seq.int((length(p)-1)/2*-1,(length(p)-1)/2)
pdf("%s_model.pdf",height=6,width=6)
plot(x,p,type="l",col=c("red"),main="Peak Model",xlab="Distance to the middle",ylab="Percentage")
lines(x,m,col=c("blue"))
lines(x,s,col=c("black"))
abline(v=median(x[p==max(p)]),lty=2,col=c("red"))
abline(v=median(x[m==max(m)]),lty=2,col=c("blue"))
abline(v=median(x[s==max(s)]),lty=2,col=c("black"))
legend("topleft",c("forward tags","reverse tags","shifted tags"),lty=c(1,1,1),col=c("red","blue","black"))
legend("right",c("d=%d"),bty="n")
dev.off()
""" % (','.join(map(str,norm_p)),','.join(map(str,norm_m)),','.join(map(str,norm_s)),name,d))
    rfhd.close()

def diag_write (filename, diag_result):
    ofhd_diag = open(filename,"w")
    a = diag_result[0]
    l = len(a)-2
    s = [90-x*10 for x in range(l)]
    ofhd_diag.write("FC range\t# of Peaks\tcover by sampling %s\n" % ("%\t".join (map(str,s))+"%"))
    format = "%s\t%d"+"\t%.2f"*l+"\n"
    ofhd_diag.write( "".join( [format % tuple(x) for x in diag_result])  )
    ofhd_diag.close()

# ------------------------------------
# Classes
# ------------------------------------
