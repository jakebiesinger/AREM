# Time-stamp: <2010-08-01 22:58:49 Tao Liu>

"""Module Description

Copyright (c) 2008 Yong Zhang, Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).

@status:  experimental
@version: $Revision$
@author:  Yong Zhang, Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""
import os
from math import log as mathlog
from array import array
from itertools import count as itertools_count

from MACS14.OutputWriter import zwig_write
from MACS14.IO.FeatIO import PeakIO,WigTrackI,BinKeeperI
from MACS14.Prob import poisson_cdf,poisson_cdf_inv
from MACS14.Constants import *

class PeakDetect:
    """Class to do the peak calling.

    e.g:
    >>> from MACS14.PeakDetect import PeakDetect
    >>> pd = PeakDetect(treat=treatdata, control=controldata, pvalue=pvalue_cutoff, d=100, scan_window=200, gsize=3000000000)
    >>> pd.call_peaks()
    >>> print pd.toxls()
    """
    def __init__ (self,opt=None,treat=None, control=None):
        """Initialize the PeakDetect object.

        """
        self.opt  = opt
        self.info = opt.info
        self.debug = opt.debug
        self.warn = opt.warn

        self.treat = treat
        self.control = control
        self.ratio_treat2control = None
        self.peaks = None
        self.final_peaks = None
        self.final_negative_peaks = None

        self.femax = opt.femax
        self.femin = opt.femin
        self.festep = opt.festep
                
        self.pvalue = opt.log_pvalue
        self.d = opt.d
        self.shift_size = self.d/2
        self.scan_window = opt.scanwindow
        self.gsize = opt.gsize
        
        self.nolambda = opt.nolambda

        self.sregion = opt.smalllocal
        self.lregion = opt.largelocal

        if (self.nolambda):
            self.info("#3 !!!! DYNAMIC LAMBDA IS DISABLED !!!!")
        self.diag = opt.diag
        self.save_wig = opt.store_wig
        #self.save_score = opt.store_score
        self.zwig_tr = opt.zwig_tr
        self.zwig_ctl= opt.zwig_ctl
        
        # params for AREM
        self.call_peak_extend_reads = False

    def call_peaks (self):
        """Call peaks function.

        Scan the whole genome for peaks. RESULTS WILL BE SAVED IN
        self.final_peaks and self.final_negative_peaks.
        """
        if self.control:                # w/ control
            self.peaks = self.__call_peaks_w_control ()
        else:                           # w/o control
            self.peaks = self.__call_peaks_wo_control ()
        return None

    def diag_result (self):
        """Run the diagnosis process on sequencing saturation.
        
        """
        if not self.diag:
            return None
        if self.control:                # w/ control
            return self.__diag_w_control()
        else:                           # w/o control
            return self.__diag_wo_control()

    def toxls (self):
        """Save the peak results in a tab-delimited plain text file
        with suffix .xls.
        
        """
        text = ""
        if self.control and self.peaks:
            text += "\t".join(("chr","start", "end",  "length",  "summit", "tags", "-10*log10(pvalue)", "fold_enrichment", "FDR(%)"))+"\n"
        elif self.peaks:
            text += "\t".join(("chr","start", "end",  "length",  "summit", "tags", "-10*log10(pvalue)", "fold_enrichment"))+"\n"
        else:
            return None
        
        chrs = self.peaks.keys()
        chrs.sort()
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                text += "%s\t%d\t%d\t%d" % (chrom,peak[0]+1,peak[1],peak[2])
                peak_summit_relative_pos = peak[3]-peak[0]
                text += "\t%d" % (peak_summit_relative_pos)
                text += "\t%d\t%.2f" % (peak[5],peak[6])
                text += "\t%.2f" % (peak[7])
                if self.control:
                    if peak[8]>=100:
                        text += "\t100"
                    else:
                        text += "\t%.2f" % (peak[8])
                text+= "\n"
        return text

    def neg_toxls (self):
        text = ""
        text += "\t".join(("chr","start", "end",  "length",  "summit", "tags", "-10*log10(pvalue)","fold_enrichment"))+"\n"
        chrs = self.final_negative_peaks.keys()
        chrs.sort()
        for chrom in chrs:
            for peak in self.final_negative_peaks[chrom]:
                text += "%s\t%d\t%d\t%d" % (chrom,peak[0]+1,peak[1],peak[2])
                peak_summit_relative_pos = peak[3]-peak[0]+1
                text += "\t%d" % (peak_summit_relative_pos)
                text += "\t%d\t%.2f" % (peak[5],peak[6])
                text += "\t%.2f" % (peak[7])
                text+= "\n"
        return text

    def tobed (self):
        text = ""
        chrs = self.peaks.keys()
        chrs.sort()
        n = 0
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                n += 1
                text+= "%s\t%d\t%d\tMACS_peak_%d\t%.2f\n" % (chrom,peak[0],peak[1],n,peak[6])
        return text

    def summitsToBED (self):
        text = ""
        chrs = self.peaks.keys()
        chrs.sort()
        n = 0
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                n += 1
                text+= "%s\t%d\t%d\tMACS_peak_%d\t%.2f\n" % (chrom,peak[3]-1,peak[3],n,peak[4])
        return text

    def __add_fdr (self, final, negative): 
        """
        A peak info type is a: dictionary

        key value: chromosome

        items: (peak start,peak end, peak length, peak summit, peak
        height, number of tags in peak region, peak pvalue, peak
        fold_enrichment, fdr) <-- tuple type
        """
        pvalue2fdr = {}
        pvalues_final = []
        pvalues_negative = []
        chrs = final.keys()
        a = pvalues_final.append
        for chrom in chrs:
            for i in final[chrom]:
                a(i[6]) # i[6] is pvalue in peak info
                pvalue2fdr[i[6]]=None
        chrs = negative.keys()
        a = pvalues_negative.append
        for chrom in chrs:
            for i in negative[chrom]:
                a(i[6])
        pvalues_final.sort(reverse=True)
        pvalues_final_l = len(pvalues_final)
        pvalues_negative.sort(reverse=True)
        pvalues_negative_l = len(pvalues_negative)        
        pvalues = pvalue2fdr.keys()
        pvalues.sort(reverse=True)
        index_p2f_pos = 0
        index_p2f_neg = 0
        for p in pvalues:
            while index_p2f_pos<pvalues_final_l and p<=pvalues_final[index_p2f_pos]:
                index_p2f_pos += 1
            n_final = index_p2f_pos

            while  index_p2f_neg<pvalues_negative_l and p<=pvalues_negative[index_p2f_neg]:
                index_p2f_neg += 1
            n_negative = index_p2f_neg
            pvalue2fdr[p] = 100.0 * n_negative / n_final

        new_info = {}
        chrs = final.keys()
        for chrom in chrs:
            new_info[chrom] = []
            for i in final[chrom]:
                tmp = list(i)
                tmp.append(pvalue2fdr[i[6]])
                new_info[chrom].append(tuple(tmp))      # i[6] is pvalue in peak info
        return new_info

    def __call_peaks_w_control (self):
        """To call peaks with control data.

        A peak info type is a: dictionary

        key value: chromosome

        items: (peak start,peak end, peak length, peak summit, peak
        height, number of tags in peak region, peak pvalue, peak
        fold_enrichment) <-- tuple type
        """
        self.lambda_bg = float(self.scan_window)*self.treat.total/self.gsize
        self.debug("#3 background lambda: %.2f " % (self.lambda_bg))
        self.min_tags = poisson_cdf_inv(1-pow(10,self.pvalue/-10),self.lambda_bg)+1
        self.debug("#3 min tags: %d" % (self.min_tags))

        self.ratio_treat2control = float(self.treat.total)/self.control.total
        if self.ratio_treat2control > 2 or self.ratio_treat2control < 0.5:
            self.warn("Treatment tags and Control tags are uneven! FDR may be wrong!")
        self.info("#3 shift treatment data")
        self.__shift_trackI(self.treat)
        self.info("#3 merge +/- strand of treatment data")

        self.treat.merge_plus_minus_locations_naive ()

        self.debug("#3 after shift and merging, tags: %d" % (self.treat.total))
        if self.save_wig:
            self.info("#3 save the shifted and merged tag counts into wiggle file...")
            #build wigtrack
            #if self.save_wig:
            #    treatwig = self.__build_wigtrackI(self.treat,space=self.opt.space)
            if self.opt.wigextend:
                zwig_write(self.treat,self.opt.wig_dir_tr,self.zwig_tr,self.opt.wigextend,log=self.info,space=self.opt.space,single=self.opt.singlewig)
            else:
                zwig_write(self.treat,self.opt.wig_dir_tr,self.zwig_tr,self.d,log=self.info,space=self.opt.space,single=self.opt.singlewig)
        self.info("#3 call peak candidates")
        peak_candidates = self.__call_peaks_from_trackI (self.treat, self.treat.prob_aligns)
        
        self.info("#3 shift control data")
        self.info("#3 merge +/- strand of control data")
        self.__shift_trackI(self.control)
        self.control.merge_plus_minus_locations_naive ()

        self.debug("#3 after shift and merging, tags: %d" % (self.control.total))
        if self.save_wig:
            self.info("#3 save the shifted and merged tag counts into wiggle file...")
            #build wigtrack
            #if self.save_score:
            #    controlbkI = self.__build_binKeeperI(self.control,space=self.opt.space)
            if self.opt.wigextend:
                zwig_write(self.control,self.opt.wig_dir_ctl,self.zwig_ctl,self.opt.wigextend,log=self.info,space=self.opt.space,single=self.opt.singlewig)
            else:
                zwig_write(self.control,self.opt.wig_dir_ctl,self.zwig_ctl,self.d,log=self.info,space=self.opt.space,single=self.opt.singlewig)
        self.info("#3 call negative peak candidates")
        negative_peak_candidates = self.__call_peaks_from_trackI (self.control, self.control.prob_aligns)
        
        if self.treat.total_multi > 0 and not self.opt.no_EM:
            self.info("#3.5 Perform EM on treatment multi reads")
            self.__align_by_EM(self.treat, self.control, peak_candidates, self.ratio_treat2control, fake_when_missing=True)

        # build score
        #if self.save_score:
        #    self.info("#4 build scores")
        #    scoreswig = self.__build_score_wigtrackI(treatwig,controlbkI,self.d,space=self.opt.space,bglambda=self.lambda_bg)
        #    zwigfile = file(self.opt.name+".score.wig","w")
        #    self.info("#4.1 save scores to wiggle file")            
        #    scoreswig.write_wig(zwigfile,"score")
        #    self.info("compress the wiggle file using gzip...")
        #    os.system("gzip "+self.opt.name+".score.wig")
        
        self.info("#3 use control data to filter peak candidates...")
        self.final_peaks = self.__filter_w_control(peak_candidates,self.treat,self.control, self.ratio_treat2control,fake_when_missing=True)
        self.info("#3 find negative peaks by swapping treat and control")

        self.final_negative_peaks = self.__filter_w_control(negative_peak_candidates,self.control,self.treat, 1.0/self.ratio_treat2control,fake_when_missing=True)
        return self.__add_fdr (self.final_peaks, self.final_negative_peaks)

    def __call_peaks_wo_control (self):
        """To call peaks w/o control data.

        """
        self.lambda_bg = float(self.scan_window)*self.treat.total/self.gsize
        self.debug("#3 background lambda: %.2f " % (self.lambda_bg))
        self.min_tags = poisson_cdf_inv(1-pow(10,self.pvalue/-10),self.lambda_bg)+1
        self.debug("#3 min tags: %d" % (self.min_tags))

        self.info("#3 shift treatment data")
        self.__shift_trackI(self.treat)
        self.info("#3 merge +/- strand of treatment data")
        self.treat.merge_plus_minus_locations_naive ()

        self.debug("#3 after shift and merging, tags: %d" % (self.treat.total))
        if self.save_wig:
            self.info("#3 save the shifted and merged tag counts into wiggle file...")
            if self.opt.wigextend:
                zwig_write(self.treat,self.opt.wig_dir_tr,self.zwig_tr,self.opt.wigextend,log=self.info,space=self.opt.space,single=self.opt.singlewig)
            else:
                zwig_write(self.treat,self.opt.wig_dir_tr,self.zwig_tr,self.d,log=self.info,space=self.opt.space,single=self.opt.singlewig)
        self.info("#3 call peak candidates")
        peak_candidates = self.__call_peaks_from_trackI (self.treat, self.treat.prob_aligns)
        self.info("#3 use self to calculate local lambda and  filter peak candidates...")
        self.final_peaks = self.__filter_w_control(peak_candidates,self.treat,self.treat,1,pass_sregion=True)
        return self.final_peaks

    def __print_peak_info (self, peak_info):
        """Print out peak information.

        A peak info type is a: dictionary

        key value: chromosome

        items: (peak start,peak end, peak length, peak summit, peak
        height, number of tags in peak region, peak pvalue, peak
        fold_enrichment) <-- tuple type
        """
        chrs = peak_info.keys()
        chrs.sort()
        for chrom in chrs:
            peak_list = peak_info[chrom]
            for peak in peak_list:
                print ( chrom+"\t"+"\t".join(map(str,peak)) )

    def __filter_w_control (self, peak_info, treatment, control, treat2control_ratio, pass_sregion=False, write2wig= False, fake_when_missing=False ):
        """Use control data to calculate several lambda values around
        1k, 5k and 10k region around peak summit. Choose the highest
        one as local lambda, then calculate p-value in poisson
        distribution.
        
        New: if clip_extra_reads, we will try to improve the p-value by clipping
        reads from both ends of the peak. This will get rid of low-probablility
        multi reads.

        Return value type in this format:
        a dictionary
        key value : chromosome
        items : array of (peak_start,peak_end,peak_length,peak_summit,peak_height,peak_num_tags,peak_pvalue,peak_fold_enrichment)
        
        """
        final_peak_info = {}
        chrs = peak_info.keys()
        chrs.sort()
        total = 0
        get_t_read_prob = lambda read: 1. if type(read) is not tuple else treatment.prob_aligns[read[1]]
        get_c_read_prob = lambda read: 1. if type(read) is not tuple else control.prob_aligns[read[1]]
        for chrom in chrs:
            self.debug("#3 Chromosome %s" % (chrom))
            n_chrom = 0
            final_peak_info[chrom] = []
            peak_list = peak_info[chrom]
            try:
                (ctags,tmp) = control.get_locations_by_chr(chrom)
            except:
                self.warn("Missing %s data, skip it..." % (chrom))
                if fake_when_missing:
                    ctags = [-1,]
                    self.warn("Fake a tag at %s:%d" % (chrom,-1))
                    tmp=[]
                else:
                    continue
            try:
                (ttags,tmp) = treatment.get_locations_by_chr(chrom)
            except:
                self.warn("Missing %s data, skip it..." % (chrom))
                if fake_when_missing:
                    ttags = [-1,]
                    self.warn("Fake a tag at %s:%d" % (chrom,-1))
                    tmp=[]
                else:
                    continue
                
            index_ctag = 0      # index for control tags
            index_ttag = 0      # index for treatment tags
            flag_find_ctag_locally = False
            flag_find_ttag_locally = False            
            prev_index_ctag = 0
            prev_index_ttag = 0            
            len_ctags =len(ctags)
            len_ttags =len(ttags)
            print '# candidates:', len(peak_list)
            for i in range(len(peak_list)):
                #(peak_start,peak_end,peak_length,peak_summit,peak_height,peak_num_tags) = peak_list[i]
                #(peak_start,peak_end,peak_length,peak_summit,peak_height,peak_num_tags, peak_indices) = peak_list[i]
                #peak_start,peak_end,peak_length,peak_summit,peak_height,number_cpr_tags, peak_num_tags, cpr_indices = peak_list[i]
                (peak_start,peak_end,peak_length,peak_summit,peak_height, peak_num_tags) = peak_list[i]

                #window_size_4_lambda = min(self.first_lambda_region,max(peak_length,self.scan_window))
                window_size_4_lambda = max(peak_length,self.scan_window)
                lambda_bg = self.lambda_bg/self.scan_window*window_size_4_lambda                
                if self.nolambda:
                    # skip local lambda
                    print 'skipping local lambda'
                    local_lambda = lambda_bg
                    tlambda_peak = float(peak_num_tags)/peak_length*window_size_4_lambda
                else:
                    left_peak = peak_start+self.shift_size # go to middle point of the first fragment
                    right_peak = peak_end-self.shift_size  # go to middle point of the last fragment
                    left_lregion = peak_summit-self.lregion/2
                    left_sregion = peak_summit-self.sregion/2
                    right_lregion = peak_summit+self.lregion/2
                    right_sregion = peak_summit+self.sregion/2
                    #(cnum_10k,cnum_5k,cnum_1k,cnum_peak) = (0,0,0,0)
                    #(tnum_10k,tnum_5k,tnum_1k,tnum_peak) = (0,0,0,0)
                    (cnum_sregion, cnum_lregion, cnum_peak, tnum_sregion, tnum_lregion, tnum_peak) = (0,0,0,0,0,0)
                    cnum_peak_total, tnum_peak_total = 0,0
                    #smallest = min(left_peak,left_10k,left_5k,left_1k)
                    #largest = max(right_peak,right_10k,right_5k,right_1k)

                    #print 'index_ctag: %s, ctags[j] %s' % (index_ctag, ctags[index_ctag])
                    while index_ctag < len_ctags:
                        if get_read_start(ctags[index_ctag]) < left_lregion:
                            # go to next control tag
                            index_ctag+=1
                        elif right_lregion < get_read_start(ctags[index_ctag]) \
                            or index_ctag + 1 >= len_ctags:
                            # finalize and go to next peak region
                            flag_find_ctag_locally = False
                            index_ctag = prev_index_ctag 
                            break
                        else:
                            if not flag_find_ctag_locally:
                                flag_find_ctag_locally = True
                                prev_index_ctag = index_ctag
                            p = get_read_start(ctags[index_ctag])
                            prob = get_c_read_prob(ctags[index_ctag])
                            if left_peak <= p <= right_peak:
                            #if peak_start <= p <= peak_end:  # Jake-- wouldn't this make more sense?
                                cnum_peak += prob
                                cnum_peak_total += 1
                            if left_sregion <= p <= right_sregion:
                                cnum_sregion += prob
                                cnum_lregion += prob
                            else:
                                cnum_lregion += prob
                            index_ctag += 1 # go to next tag

                    inds_in_peak = []
                    while index_ttag < len_ttags:
                        if get_read_start(ttags[index_ttag]) < left_lregion:
                            # go to next treatment tag
                            index_ttag+=1
                        elif right_lregion < get_read_start(ttags[index_ttag]) \
                            or index_ttag + 1 >= len_ttags:
                            # finalize and go to next peak region
                            flag_find_ttag_locally = False
                            index_ttag = prev_index_ttag
                            break
                        else:
                            if not flag_find_ttag_locally:
                                flag_find_ttag_locally = True
                                prev_index_ttag = index_ttag
                            p = get_read_start(ttags[index_ttag])
                            prob = get_t_read_prob(ttags[index_ttag])
                            if left_peak <= p <= right_peak:
                            #if peak_start <= p <= peak_end:  # Jake-- again, seems to be more accurate...
                                tnum_peak += prob
                                tnum_peak_total += 1
                                inds_in_peak.append(index_ttag)
                            if left_sregion <= p <= right_sregion:
                                tnum_sregion += prob
                                tnum_lregion += prob
                            else:
                                tnum_lregion += prob
                            index_ttag += 1 # go to next tag
                    clambda_peak = float(cnum_peak)/peak_length*treat2control_ratio*window_size_4_lambda
                    #clambda_10k = float(cnum_10k)/self.third_lambda_region*treat2control_ratio*window_size_4_lambda
                    clambda_lregion = float(cnum_lregion)/self.lregion*treat2control_ratio*window_size_4_lambda
                    clambda_sregion = float(cnum_sregion)/self.sregion*treat2control_ratio*window_size_4_lambda
                    tlambda_peak = float(tnum_peak)/peak_length*window_size_4_lambda
                    #tlambda_10k = float(tnum_10k)/self.third_lambda_region*window_size_4_lambda
                    tlambda_lregion = float(tnum_lregion)/self.lregion*window_size_4_lambda
                    tlambda_sregion = float(tnum_sregion)/self.sregion*window_size_4_lambda
                    #print clambda_peak, clambda_lregion, clambda_sregion, tlambda_peak, tlambda_sregion, tlambda_lregion

                    if pass_sregion:
                        # for experiment w/o control, peak region lambda and sregion region lambda are ignored!
                        local_lambda = max(lambda_bg,tlambda_lregion)
                    else:
                        # for experiment w/ control
                        local_lambda = max(lambda_bg,clambda_peak,clambda_lregion,clambda_sregion)

                p_tmp = poisson_cdf(tlambda_peak,local_lambda,lower=False)
                if p_tmp <= 0:
                    peak_pvalue = 3100
                else:
                    peak_pvalue = mathlog(p_tmp,10) * -10
                
                #print 'counts:', cnum_peak, cnum_sregion, cnum_lregion
                #print 'lambdas:', clambda_peak, clambda_sregion, clambda_lregion, lambda_bg
                
                if self.call_peak_extend_reads:
                    # call subpeaks from peak heights. Extend each read by scan_window 
                    # and require a min height corresponding to the bg_lambda.
                    # bases out from the tag center.  Then call subpeaks by looking
                    # for regions that are at least as high as the pvalue minimum.
                    # separate regions that don't pass this criteria
                    peak_tags = [ttags[j] for j in inds_in_peak]
                    #peak_height_lambda = self.lambda_bg / self.scan_window * 
                    thresh_upper = poisson_cdf_inv(1-pow(10,self.pvalue/-10),local_lambda)
                    pval_low = poisson_cdf(thresh_upper, local_lambda, lower=False)
                    pval_high = poisson_cdf(thresh_upper - 1, local_lambda, lower=False)
                    pval_mid = pow(10, self.pvalue/-10)
                    thresh_subpeak = thresh_upper + -1./(pval_high - pval_low) * (pval_mid - pval_high)  # rise of -1, run of diff(pvals)
                    #print self.min_tags, self.scan_window
                    #print peak_tags
                    #print 'min threshold for local_lambda %s, peak_length %s is %s' % (local_lambda, peak_length, thresh_subpeak)
                    subpeaks = self.__tags_call_peak_w_subpeaks(peak_tags,
                                        treatment.prob_aligns, thresh_subpeak)
                    #print subpeaks
                    for sub_p in subpeaks:
                        n_chrom += 1
                        total += 1
                        (peak_start,peak_end,peak_length,peak_summit,
                            peak_height, starts_in_peak) = sub_p
                        peak_num_tags = len(starts_in_peak)
                        peak_mass = sum([get_t_read_prob(p) for p in starts_in_peak])
                        p_tmp = poisson_cdf(peak_mass,local_lambda,lower=False)
                        if p_tmp <= 0:
                            peak_pvalue = 3100
                        else:
                            peak_pvalue = mathlog(p_tmp,10) * -10
                        peak_fold_enrichment = float(peak_height)/local_lambda*window_size_4_lambda/self.d
                        final_peak_info[chrom].append((peak_start,peak_end,peak_length,peak_summit,peak_height,peak_mass,peak_pvalue,peak_fold_enrichment))

                elif not self.opt.no_greedy_caller:
                    # build up sub peaks from the candidate we are iterating over
                    # by greedily adding tags to the current subpeak.  To avoid
                    # local minima, we always look at least scanwindow bases away
                    # get reads within first scanwindow
                    cpr_tags = [ttags[inds_in_peak[0]]]
                    cpr_mass = get_t_read_prob(ttags[inds_in_peak[0]])
                    j = 1
                    while j < len(inds_in_peak):
                        cur_dist = get_read_start(ttags[inds_in_peak[j]]) - get_read_start(cpr_tags[0])
                        if cur_dist <= self.scan_window:
                            cpr_tags.append(ttags[inds_in_peak[j]])
                            cpr_mass += get_t_read_prob(cpr_tags[-1])
                            j += 1
                        else:
                            break
                    # add reads to current peak if they improve the enrichment
                    cpr_width = get_read_start(cpr_tags[-1]) - get_read_start(cpr_tags[0])
                    lambda_width = max(cpr_width, self.scan_window)
                    cpr_pval = poisson_cdf(cpr_mass,local_lambda * lambda_width / window_size_4_lambda,lower=False)
                    middle_mass = 0
                    #print 'on first', cpr_width, cpr_mass, cpr_pval, cpr_tags
                    #print j, len(inds_in_peak)
                    while j < len(inds_in_peak):
                        #print cpr_width, cpr_mass, cpr_pval, cpr_tags
                        test_posn = get_read_start(ttags[inds_in_peak[j]])
                        cur_dist = test_posn - get_read_start(cpr_tags[-1])
                        if cur_dist <= self.scan_window:
                            test_width = test_posn - get_read_start(cpr_tags[0])
                            lambda_width = max(test_width, self.scan_window)
                            test_mass = cpr_mass + middle_mass + get_t_read_prob(ttags[inds_in_peak[j]])
                            test_pval = poisson_cdf(test_mass,local_lambda * lambda_width / window_size_4_lambda,lower=False)
                            #print 'vs.', test_width, test_mass, test_pval
                            if test_pval < cpr_pval:
                                # enrichment improved-- add the tag
                                cpr_tags.append(ttags[inds_in_peak[j]])
                                cpr_mass = test_mass
                                cpr_width = test_width
                                cpr_pval = test_pval
                                middle_mass = 0
                                #print 'accepted'
                            else:
                                middle_mass += get_t_read_prob(ttags[inds_in_peak[j]])
                                #print 'denied'
                            j += 1
                        else:
                            # call previous region as a peak
                            if cpr_pval <= 0:
                                cpr_pval = 3100
                            else:
                                cpr_pval = mathlog(cpr_pval,10) * -10
                            #print 'calling peak!', cpr_pval, '>', self.pvalue, ' ?'
                            if cpr_pval > self.pvalue:
                                p_start, p_end, p_length, p_summit, p_height = self.__tags_call_peak(cpr_tags, treatment.prob_aligns)
                                cpr_enrich = p_height / local_lambda * window_size_4_lambda / self.d
                                final_peak_info[chrom].append((p_start, p_end, p_length, p_summit, p_height, cpr_mass, cpr_pval, cpr_enrich))
                                n_chrom += 1
                                total += 1

                            # reset cpr
                            cpr_tags = [ttags[inds_in_peak[j]]]
                            cpr_mass = get_t_read_prob(ttags[inds_in_peak[j]])
                            while j < len(inds_in_peak):
                                cur_dist = get_read_start(ttags[inds_in_peak[j]]) - get_read_start(cpr_tags[0])
                                if cur_dist <= self.scan_window:
                                    cpr_tags.append(ttags[inds_in_peak[j]])
                                    cpr_mass += get_t_read_prob(cpr_tags[-1])
                                    j += 1
                                else:
                                    break
                            cpr_width = get_read_start(cpr_tags[-1]) - get_read_start(cpr_tags[0])
                            if cpr_width > 0:
                                lambda_width = max(cpr_width, self.scan_window)
                                cpr_pval = poisson_cdf(cpr_mass,local_lambda * lambda_width / window_size_4_lambda,lower=False)
                            else:
                                cpr_pval = 1.
                            middle_mass = 0
                        # clip last tag from the left if it improves enrichment
                        if len(cpr_tags) > 3:
                            test_width = get_read_start(cpr_tags[-1]) - get_read_start(cpr_tags[1])
                            test_mass = cpr_mass - get_t_read_prob(cpr_tags[0])
                            lambda_width = max(test_width, self.scan_window)
                            test_pval = poisson_cdf(test_mass,local_lambda * lambda_width / window_size_4_lambda,lower=False)
                            if test_pval < cpr_pval:
                                # enrichment improved-- add the tag
                                cpr_tags.pop(0)
                                cpr_mass = test_mass
                                cpr_width = test_width
                                cpr_pval = test_pval
                        
                    if len(cpr_tags) > 1:  # call last reads
                        if cpr_pval <= 0:
                            cpr_pval = 3100
                        else:
                            cpr_pval = mathlog(cpr_pval,10) * -10
                        #print 'calling peak!', cpr_pval, '>', self.pvalue, ' ?'
                        if cpr_pval > self.pvalue:
                            p_start, p_end, p_length, p_summit, p_height = self.__tags_call_peak(cpr_tags, treatment.prob_aligns)
                            cpr_enrich = p_height / local_lambda * window_size_4_lambda / self.d
                            final_peak_info[chrom].append((p_start, p_end, p_length, p_summit, p_height, cpr_mass, cpr_pval, cpr_enrich))
                            n_chrom += 1
                            total += 1

                elif peak_pvalue > self.pvalue:
                    n_chrom += 1
                    total += 1
                    peak_fold_enrichment = float(peak_height)/local_lambda*window_size_4_lambda/self.d
                    final_peak_info[chrom].append((peak_start,peak_end,peak_length,peak_summit,peak_height,peak_num_tags,peak_pvalue,peak_fold_enrichment))
                #else:
                #    self.debug("Reject the peak at %s:%d-%d with local_lambda: %.2f and -log10pvalue: %.2f" % (chrom,peak_start,peak_end,local_lambda,peak_pvalue))

            self.debug("#3 peaks whose pvalue < cutoff: %d" % (n_chrom))
        self.info("#3 Finally, %d peaks are called!" % (total))
        return final_peak_info

    def __call_peaks_from_trackI (self, trackI, prob_aligns):
        """Call peak candidates from trackI data. Using every tag as
        step and scan the self.scan_window region around the tag. If
        tag number is greater than self.min_tags, then the position is
        recorded.

        Return: data in this format. (peak_start,peak_end,peak_length,peak_summit,peak_height,peak_num_tags)
        """
        peak_candidates = {}
        self.debug("#3 search peak condidates...")
        chrs = trackI.get_chr_names()
        total = 0
        get_read_prob = lambda read: 1. if type(read) is not tuple else prob_aligns[read[1]]
        for chrom in chrs:
            self.debug("#3 Chromosome %s" % (chrom))
            n_chrom = 0
            peak_candidates[chrom] = []
            (tags,tmp) = trackI.get_locations_by_chr(chrom)
            len_t = len(tags)
            cpr_tags = []       # Candidate Peak Region tags
            cpr_tags.extend(tags[:self.min_tags-1])
            number_cpr_tags = self.min_tags-1
            p = self.min_tags-1 # Next Tag Index
            cpr_probs = [get_read_prob(tags[i]) for i in range(p+1)]
            #cpr_indices = range(p+1)
            while p < len_t:
                if number_cpr_tags >= self.min_tags:
                    if get_read_start(tags[p]) - get_read_start(cpr_tags[-1*self.min_tags+1]) <= self.scan_window:
                        # add next tag, if the new tag is less than self.scan_window away from previous no. self.min_tags tag
                        cpr_tags.append(tags[p])
                        cpr_probs.append(get_read_prob(tags[p]))
                        #cpr_indices.append(p)
                        number_cpr_tags += 1
                        p+=1
                    else:
                        # candidate peak region is ready, call peak...
                        (peak_start,peak_end,peak_length,peak_summit,peak_height) = self.__tags_call_peak (cpr_tags, prob_aligns)
                        peak_prob = sum(cpr_probs)
                        #peak_candidates[chrom].append((peak_start,peak_end,peak_length,peak_summit,peak_height,number_cpr_tags, peak_prob, cpr_indices))
                        peak_candidates[chrom].append((peak_start,peak_end,peak_length,peak_summit,peak_height,peak_prob))
                        cpr_tags = [tags[p]] # reset
                        cpr_probs = [get_read_prob(tags[p])]
                        #cpr_indices = [p]
                        number_cpr_tags = 1
                        total += 1
                        n_chrom += 1
                        p += 1
                else:
                    # add next tag, but if the first one in cpr_tags
                    # is more than self.scan_window away from this new
                    # tag, remove the first one from the list
                    if get_read_start(tags[p]) - get_read_start(cpr_tags[0]) >= self.scan_window:
                        cpr_tags.pop(0)
                        number_cpr_tags -= 1
                        cpr_probs.pop(0)
                        #cpr_indices.pop(0)
                    cpr_tags.append(tags[p])
                    cpr_probs.append(get_read_prob(tags[p]))
                    #cpr_indices.append(p)
                    number_cpr_tags += 1
                    p+=1
            self.debug("#3 peak candidates: %d" % (n_chrom))
        self.debug("#3 Total number of candidates: %d" % (total))
        return self.__remove_overlapping_peaks(peak_candidates)
                
    def __tags_call_peak (self, tags, prob_aligns ):
        """Project tags to a line. Then find the highest point.

        """
        start = get_read_start(tags[0])-self.d/2
        end = get_read_start(tags[-1])+self.d/2       # +1 or not?
        region_length = end - start
        line= [0]*region_length
        for tag in tags:
            t_start = get_read_start(tag)
            t_prob = 1. if type(tag) is not tuple else prob_aligns[tag[1]]
            tag_projected_start = t_start-start-self.d/2
            tag_projected_end = t_start-start+self.d/2
            for i in range(tag_projected_start,tag_projected_end):
                line[i] += t_prob
        tops = []
        top_height = 0
        for i in range(len(line)):
            if line[i] > top_height:
                top_height = line[i]
                tops = [i]
            elif line[i] == top_height:
                tops.append(i)
        peak_summit = tops[len(tops)/2]+start
        return (start,end,region_length,peak_summit,top_height)

    def __tags_call_peak_w_subpeaks (self, tags, prob_aligns, subpeak_threshold):
        """Project tags to a line, extending them by self.window_size.
        break the peak region into subpeaks by looking for non-contiguous
        regions that pass the p-value threshold.
        
        return a list of subpeaks of the format:
        [(sub_start, sub_end, sub_length, peak_summit, top_height, starts_in_peak), ...]
        """
        start = get_read_start(tags[0])-self.scan_window
        end = get_read_start(tags[-1])+self.scan_window      # +1 or not?
        region_length = end - start
        line= [0.]*region_length
        for tag in tags:
            t_start = get_read_start(tag)
            t_prob = 1. if type(tag) is not tuple else prob_aligns[tag[1]]
            tag_projected_start = t_start-start-self.scan_window
            tag_projected_end = t_start-start+self.scan_window
            for i in range(tag_projected_start,tag_projected_end):
                line[i] += t_prob
        tops = []
        top_height = 0
        in_peak = False
        subpeaks = []
        sub_start, sub_end = 0,0
        for i in range(len(line)):
            if line[i] >= subpeak_threshold:
                if line[i] > top_height:
                    top_height = line[i]
                    tops = [i]
                elif line[i] == top_height:
                    tops.append(i)
                if in_peak:  # extend the peak region
                    sub_end = start + i + 1
                else:  # start the peak region
                    sub_start, sub_end = start + i, start + i+1
                    in_peak = True
                if in_peak and i + 1 >= len(line):  # need to call the last peak
                    peak_summit = tops[len(tops)/2]+start
                    sub_length = sub_end - sub_start
                    subpeaks.append((sub_start, sub_end, sub_length, peak_summit, top_height))
                    in_peak = False
                    tops = []
                    top_height = 0
            elif in_peak:  # close the peak region
                peak_summit = tops[len(tops)/2]+start
                sub_length = sub_end - sub_start
                subpeaks.append((sub_start, sub_end, sub_length, peak_summit, top_height))
                in_peak = False
                tops = []
                top_height = 0
        # figure out which tags are in which subpeaks
        tag_i = 0
        for sub_i in range(len(subpeaks)):
            sub_start, sub_end = subpeaks[sub_i][:2]
            starts_in_peak = []
            while tag_i < len(tags) and get_read_start(tags[tag_i]) <= sub_end:
                if get_read_start(tags[tag_i]) >= sub_start:
                    starts_in_peak.append(tags[tag_i])
                tag_i += 1
            assert len(starts_in_peak) > 0
            subpeaks[sub_i] = subpeaks[sub_i] + (starts_in_peak,)
        return subpeaks

    def __shift_trackI (self, trackI):
        """Shift trackI data to right (for plus strand) or left (for
        minus strand).

        """
        chrs = trackI.get_chr_names()
        for chrom in chrs:
            tags = trackI.get_locations_by_chr(chrom)
            # plus
            for i in range(len(tags[0])):
                t = tags[0][i]
                #print t
                if type(t) is int:
                    tags[0][i] += self.shift_size
                else:
                    tags[0][i] = (t[0] + self.shift_size, t[1])
            # minus
            for i in range(len(tags[1])):
                t = tags[1][i]
                if type(t) is int:
                    tags[1][i] -= self.shift_size
                else:
                    tags[1][i] = (t[0] - self.shift_size, t[1])
        return
    
    def __build_wigtrackI (self, trackI, space=10):
        """Shift trackI then build a wigTrackI object.
        
        """
        chrs = trackI.get_chr_names()
        wigtrack = WigTrackI()
        wigtrack.span = space
        d = self.d
        step = 10000000 + 2*d
        
        for chrom in chrs:
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
                    # keep this tag for next window
                    for i in xrange(d,step-d,space):
                        if window_counts[i] == 0:
                            pass
                        else:
                            wigtrack.add_loc(chrom,i+startp+1,window_counts[i])
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
                    wigtrack.add_loc(chrom,i+startp+1,window_counts[i])                    
        return wigtrack

    def __diag_w_control (self):
        # sample
        sample_peaks = {}
        for i in xrange(90,10,-10):
            self.info("#3 diag: sample %d%%" % i)
            sample_peaks[i]=self.__diag_peakfinding_w_control_sample(float(i)/(i+10))
        return self.__overlap (self.final_peaks, sample_peaks,top=90,bottom=10,step=-10)

    def __diag_peakfinding_w_control_sample (self, percent):
        self.treat.sample(percent) # because sampling is after
                                   # shifting, track.total is used
                                   # now.
        self.control.sample(percent)
        ratio_treat2control = float(self.treat.total)/self.control.total

        self.lambda_bg = float(self.scan_window)*self.treat.total/self.gsize # bug fixed...
        self.min_tags = poisson_cdf_inv(1-pow(10,self.pvalue/-10),self.lambda_bg)+1

        self.debug("#3 diag: after shift and merging, treat: %d, control: %d" % (self.treat.total,self.control.total))
        self.info("#3 diag: call peak candidates")
        peak_candidates = self.__call_peaks_from_trackI (self.treat, self.treat.prob_aligns)

        self.info("#3 diag: call negative peak candidates")
        negative_peak_candidates = self.__call_peaks_from_trackI (self.control, self.control.prob_aligns)
        
        self.info("#3 diag: use control data to filter peak candidates...")
        final_peaks_percent = self.__filter_w_control(peak_candidates,self.treat,self.control, ratio_treat2control)
        return final_peaks_percent
        
    def __diag_wo_control (self):
        # sample
        sample_peaks = {}
        for i in xrange(90,10,-10):
            self.info("#3 diag: sample %d%%" % i)
            sample_peaks[i]=self.__diag_peakfinding_wo_control_sample(float(i)/(i+10))
        return self.__overlap (self.final_peaks, sample_peaks,top=90,bottom=10,step=-10)

    def __diag_peakfinding_wo_control_sample (self, percent):

        self.lambda_bg = float(self.scan_window)*self.treat.total/self.gsize # bug fixed...
        self.min_tags = poisson_cdf_inv(1-pow(10,self.pvalue/-10),self.lambda_bg)+1

        self.treat.sample(percent)
        self.debug("#3 diag: after shift and merging, tags: %d" % (self.treat.total))
        self.info("#3 diag: call peak candidates")
        peak_candidates = self.__call_peaks_from_trackI (self.treat, self.treat.prob_aligns)
        self.info("#3 diag: use self to calculate local lambda and  filter peak candidates...")
        final_peaks_percent = self.__filter_w_control(peak_candidates,self.treat,self.treat,1,pass_sregion=True) # bug fixed...
        return final_peaks_percent

    def __overlap (self, gold_peaks, sample_peaks, top=90,bottom=10,step=-10):
        """Calculate the overlap between several fe range for the
        golden peaks set and results from sampled data.
        
        """
        gp = PeakIO()
        gp.init_from_dict(gold_peaks)
        if self.femax:
            femax = min(self.femax, (int(gp.max_fold_enrichment())//self.festep+1)*self.festep)
        else:
            femax = (int(gp.max_fold_enrichment())//self.festep+1)*self.festep
        femin = self.femin
        diag_result = []
        for f in xrange(femin, femax, self.festep):
            
            fe_low = f
            fe_up = f + self.festep
            self.debug("#3 diag: fe range = %d -- %d" % (fe_low, fe_up))
            
            r = self.__overlap_fe(gold_peaks, sample_peaks, fe_low, fe_up, top, bottom, step)
            if r:
                diag_result.append(r)
        return diag_result

    def __overlap_fe (self, gold_peaks, sample_peaks, fe_low, fe_up, top, bottom, step):
        ret = ["%d-%d" % (fe_low,fe_up)]
        gp = PeakIO()
        gp.init_from_dict(gold_peaks)
        gp.filter_fc(fe_low,fe_up)
        gptotal =  gp.total()
        if gptotal <= 0:
            return None

        ret.append(gptotal)
        for i in xrange(top,bottom,step):
            p = PeakIO()
            p.init_from_dict(sample_peaks[i])
            percent = 100.0*gp.overlap_with_other_peaks(p)/gptotal
            ret.append(percent)
            del p
        return ret


    def __remove_overlapping_peaks (self, peaks ):
        """peak_candidates[chrom] = [(peak_start,peak_end,peak_length,peak_summit,peak_height,number_cpr_tags)...]

        """
        new_peaks = {}
        chrs = peaks.keys()
        chrs.sort()
        for chrom in chrs:
            new_peaks[chrom]=[]
            n_append = new_peaks[chrom].append
            prev_peak = None
            peaks_chr = peaks[chrom]
            for i in xrange(len(peaks_chr)):
                if not prev_peak:
                    prev_peak = peaks_chr[i]
                    continue
                else:
                    if peaks_chr[i][0] <= prev_peak[1]:
                        s_new_peak = prev_peak[0]
                        e_new_peak = peaks_chr[i][1]
                        l_new_peak = e_new_peak-s_new_peak
                        if peaks_chr[i][4] > prev_peak[4]:
                            summit_new_peak = peaks_chr[i][3]
                            h_new_peak = peaks_chr[i][4]
                        else:
                            summit_new_peak = prev_peak[3]
                            h_new_peak = prev_peak[4]
                        prev_peak = (s_new_peak,e_new_peak,l_new_peak,summit_new_peak,h_new_peak,peaks_chr[i][5]+prev_peak[5])
                    else:
                        n_append(prev_peak)
                        prev_peak = peaks_chr[i]
            if prev_peak:
                n_append(prev_peak)
        return new_peaks


    def __align_by_EM(self, treatment, control, init_regions, treat2control_ratio, pass_sregion=False, fake_when_missing=False):
        """
        Align the multi reads in treatment using expectation-maximization.
        
        This process progressively assigns reads with multiple mappings to
        the most enriched of the possible mappings.  We require a list of
        candidate regions (init_regions) before getting started. This list is
        generated using the MACS peak caller, considering each alignment of all
        multi-reads to have as much weight as a unique read. Although these
        initial peaks will on average be larger and more diffuse than the
        finally-called peaks, we are trying to catch clustering multi-reads.
        
        """
        # get indices and lambdas for candidate regions
        all_peak_inds, all_peak_lambdas = self.__get_all_peak_lambdas(init_regions,
            treatment, control, treat2control_ratio, pass_sregion, fake_when_missing)
        # filter out non-multiread peaks and save only counts of unique reads + multi indices
        min_score = self.opt.min_score if self.opt.min_score is not None else 1e-3
        max_score = self.opt.max_score if self.opt.max_score is not None else 2e3
        final_regions = {}
        peak_posns = {}
        show_graphs = self.opt.show_graphs
        #if show_graphs:
        in_candidate = [False] * len(treatment.prob_aligns)  # if peak is in cand region
        for chrom in all_peak_inds:
            try:
                (ttags,tmp) = treatment.get_locations_by_chr(chrom)
            except:
                self.warn("Missing %s data, skip it..." % (chrom))
                if fake_when_missing:
                    ttags = [-1,]
                    self.warn("Fake a tag at %s:%d" % (chrom,-1))
                    tmp=[]
                else:
                    continue
            for i in range(len(all_peak_inds[chrom])):
                local_lambda = all_peak_lambdas[chrom][i]
                peak_inds = all_peak_inds[chrom][i]
                unique_count = 0
                multi_inds = []
                if chrom not in peak_posns:
                    peak_posns[chrom] = []
                peak_posns[chrom].append((get_read_start(ttags[peak_inds[0]]),
                                          get_read_start(ttags[peak_inds[-1]])))
                for ind in peak_inds:
                    if type(ttags[ind]) is tuple:
                        multi_inds.append(ttags[ind][1])  # get the index into the prob array
                        #if show_graphs:
                        in_candidate[ttags[ind][1]] = True
                    else:
                        unique_count += 1
                if len(multi_inds) == 0:
                    continue  # no multi-reads in this peak.  Skip it
                else:
                    if chrom not in final_regions:
                        final_regions[chrom] = []
                    final_regions[chrom].append((local_lambda, unique_count, multi_inds))
        print 'total_multi: ', treatment.total_multi
        print 'total number of peaks in candidate regions:', sum(1 for inc in in_candidate if inc)
        
        # for each iteration
        #if False:
        if show_graphs:
            self.__plot_EM_state(0, final_regions, peak_posns, in_candidate)
        prev_entropy = None
        log_base = self.opt.enrich_log_base
        for iteration in itertools_count(1):  # until convergence
            cur_entropy = list(self.__get_read_entropy(self.treat))
            if prev_entropy is not None:  # check for convergence
                denom = sum(ent ** 2 for ent in cur_entropy)
                if denom == 0.0:
                    difference = 0.0
                else:
                    difference = sum([(cur_entropy[i] - prev_entropy[i])**2 
                        for i in xrange(len(cur_entropy))])
                    difference = difference / denom
                self.info("Entropy difference is %s" % difference)
                if difference < self.opt.min_change:
                    self.info("Convergence criteria reached after %s iterations!" % iteration)
                    break
            self.info("#3.%s iterate AREM" % iteration)
            # calculate the enrichment of each candidate peak
            for chrom in sorted(final_regions.keys()):
                for local_lambda, unique_mass, multi_inds in final_regions[chrom]:
                    multi_mass = sum(treatment.prob_aligns[i] for i in multi_inds)
                    pvalue = poisson_cdf(multi_mass + unique_mass, local_lambda,
                                         lower=False)
                    if pvalue <= 0:
                        score = max_score
                    else:
                        score = max(-mathlog(pvalue, log_base), min_score)
                        score = min(score, max_score)
                    for i in multi_inds:
                        treatment.enrich_scores[i] = score
            # normalize the alignment probabilities for all multireads
            for i in range(len(treatment.group_starts)):
                group_start = treatment.group_starts[i]
                if i < len(treatment.group_starts) - 1:
                    group_end = treatment.group_starts[i+1]
                else:
                    group_end = len(treatment.prob_aligns)
                group_range = range(group_start, group_end)
                #print [treatment.enrich_scores[j] for j in group_range]
                #print [treatment.prob_aligns[j] for j in group_range]
                enrich_total = sum(treatment.enrich_scores[j] for j in group_range)
                for j in group_range:
                    treatment.prob_aligns[j] = treatment.enrich_scores[j] / enrich_total
            if show_graphs:
                self.__plot_EM_state(iteration, final_regions, peak_posns, in_candidate)
            prev_entropy = cur_entropy
            # rinse and repeat (until convergence)
    
    def __plot_EM_state(self, iteration, final_regions, peak_posns, in_candidate, output_summary=False):
        '''Plot the current data state. These may include:
           Entropy Histogram, CDF of enrichment scores and alignment probabilities,
           ratio of FG to BG scores or probabilities
        '''
        if output_summary:
            with open(self.opt.name + '_EMpeaks_%s.txt' % iteration, 'w') as outfile:
                # final_regions[chrom].append((local_lambda, unique_count, multi_inds))
                outfile.write('\t'.join(['chrom', 'start', 'stop', 'length', 'local_lamba',
                                         'unique_count', 'total_mass']) + '\n')
                for chrom in final_regions:
                    for i, data in enumerate(final_regions[chrom]):
                        peak_mass = sum(self.treat.prob_aligns[j] for j in data[2])
                        peak_mass += data[1]
                        peak_length = peak_posns[chrom][i][1] - peak_posns[chrom][i][0]
                        line = (chrom,) + peak_posns[chrom][i] +  (peak_length, data[0], data[1], peak_mass)
                        outfile.write('\t'.join(map(str, line)) + '\n')
        self.__plot_entropy_hist(iteration, self.treat)
        self.__plot_enrichment_CDF(iteration, self.treat)
        self.__plot_probability_CDF(iteration, self.treat)
        #self.__plot_probability_mass_ratio(iteration, self.treat, self.control)
        self.__plot_max_probability(iteration, self.treat)

    def __get_read_entropy(self, treatment, normed=True, minScore=None):
        'generator for entropy in read alignments'
        for i in range(len(treatment.group_starts)):
            group_start = treatment.group_starts[i]
            if i < len(treatment.group_starts) - 1:
                group_end = treatment.group_starts[i+1]
            else:
                group_end = len(treatment.prob_aligns)
            group_range = range(group_start, group_end)
            if minScore is None:
                probs = [treatment.prob_aligns[j] for j in group_range]
            else:
                scores = [treatment.enrichScore[j] - minScore for j in group_range]
                scoreTotal = sum(scores)
                if scoreTotal < 0.:
                    scoreTotal = 0.
                # renormalize, or if no alignments had score > minScore, consider them uniform
                if scoreTotal > 1e-5:
                    probs = [score / scoreTotal for score in scores]
                else:
                    probs = [1./len(group_range)] * len(group_range)
            entropies = [p * mathlog(p) if p > 0 else 0. for p in probs]
            if normed:
                yield -sum(entropies) / mathlog(len(group_range))
            else:
                yield -sum(entropies)

    def __plot_entropy_hist(self, iteration, treatment):
        from matplotlib import pyplot
        entropy = list(self.__get_read_entropy(treatment))
        #outfile = open('entropy.%s.txt' % iteration, 'wb')
        #for value in entropy:
            #outfile.write(str(value)+'\n')
        #outfile.close()
        n, bins, patches = pyplot.hist(entropy, 50, facecolor='black', alpha=1)
        #pyplot.xticks( scipy.arange(0,1.1,0.1) )
        pyplot.xlim([0,1])
        pyplot.xlabel('Relative Entropy')
        pyplot.ylabel('Number of reads')
        pyplot.title('Multihit entropy distribution for %s at i=%s' % (
                        self.opt.name,iteration))
        pyplot.savefig(self.opt.name + '_entropy_%s.png' % iteration)
        pyplot.close()

    def __plot_enrichment_CDF(self, iteration, treatment):
        from matplotlib import pyplot
        pyplot.hist(treatment.enrich_scores, bins=100, normed=True, 
                    cumulative=True, histtype='step') 
        #outfile = open('enrichScore.%s.txt' % iteration, 'wb')
        #for value in treatment.enrich_scores:
            #outfile.write(str(value)+'\n')
        #outfile.close()
        pyplot.ylim([0,1])
        pyplot.xlabel('Enrichment Score')
        pyplot.ylabel('fraction of data')
        pyplot.title('CDF of Enrichment Scores for %s at i=%s' % (self.opt.name,
                                                                  iteration))
        pyplot.savefig(self.opt.name + '_CDF_enrichment_%s.png' % iteration)
        pyplot.close()
    
    def __plot_probability_CDF(self, iteration, treatment):
        from matplotlib import pyplot
        #outfile = open('alignProbs.%s.txt' % iteration, 'wb')
        #for value in treatment.prob_aligns:
            #outfile.write(str(value)+'\n')
        #outfile.close()
        pyplot.hist(treatment.prob_aligns, bins=100, normed=True, 
                     cumulative=True, histtype='step') 
        #pyplot.xticks( scipy.arange(0,1.1,0.1) )
        pyplot.xlim([0,1])
        pyplot.ylim([0,1])
        pyplot.xlabel('Alignment Probability')
        pyplot.ylabel('Fraction of data')
        pyplot.title('CDF of alignment probabilities for %s at i=%s' % (
                        self.opt.name,iteration))
        pyplot.savefig(self.opt.name + '_CDF_prob_%s.png' % iteration)
        pyplot.close()
    
    def __plot_max_probability(self, iteration, treatment):
        from matplotlib import pyplot
        max_probs = []
        for i in range(len(treatment.group_starts)):
            group_start = treatment.group_starts[i]
            if i < len(treatment.group_starts) - 1:
                group_end = treatment.group_starts[i+1]
            else:
                group_end = len(treatment.prob_aligns)
            group_range = range(group_start, group_end)
            max_probs.append(max([treatment.prob_aligns[i] for i in group_range]))
        pyplot.hist(max_probs, bins=50, facecolor='black', alpha=1)
        pyplot.xlim([0,1])
        pyplot.xlabel('Highest Alignment Probability')
        pyplot.ylabel('Count')
        pyplot.title('Highest read alignment probability for %s at i=%s' % (self.opt.name, iteration))
        pyplot.savefig(self.opt.name + '_max_prob_%s.png' % iteration)
        pyplot.close()


    def __get_all_peak_lambdas(self, peak_info, treatment, control, treat2control_ratio, pass_sregion=False, fake_when_missing=False):
        """
        from MACS: calculate the local (max) lambda for each peak.
        Also returns all tag indices within each peak
        
        """
        chroms = sorted(peak_info.keys())
        all_peak_inds = {}
        all_local_lambdas = {}
        get_t_read_prob = lambda read: 1. if type(read) is not tuple else treatment.prob_aligns[read[1]]
        get_c_read_prob = lambda read: 1. if type(read) is not tuple else control.prob_aligns[read[1]]

        for chrom in chroms:
            peak_list = peak_info[chrom]
            try:
                (ctags,tmp) = control.get_locations_by_chr(chrom)
            except:
                self.warn("Missing %s data, skip it..." % (chrom))
                if fake_when_missing:
                    ctags = [-1,]
                    self.warn("Fake a tag at %s:%d" % (chrom,-1))
                    tmp=[]
                else:
                    continue
            try:
                (ttags,tmp) = treatment.get_locations_by_chr(chrom)
            except:
                self.warn("Missing %s data, skip it..." % (chrom))
                if fake_when_missing:
                    ttags = [-1,]
                    self.warn("Fake a tag at %s:%d" % (chrom,-1))
                    tmp=[]
                else:
                    continue
                
            index_ctag = 0      # index for control tags
            index_ttag = 0      # index for treatment tags
            flag_find_ctag_locally = False
            flag_find_ttag_locally = False            
            prev_index_ctag = 0
            prev_index_ttag = 0            
            len_ctags =len(ctags)
            len_ttags =len(ttags)            
            for i in range(len(peak_list)):
                #(peak_start,peak_end,peak_length,peak_summit,peak_height,peak_num_tags) = peak_list[i]
                #peak_start,peak_end,peak_length,peak_summit,peak_height,number_cpr_tags, peak_num_tags, cpr_tags = peak_list[i]
                peak_start,peak_end,peak_length,peak_summit,peak_height, peak_num_tags = peak_list[i]
        
                #window_size_4_lambda = min(self.first_lambda_region,max(peak_length,self.scan_window))
                window_size_4_lambda = max(peak_length,self.scan_window)
                #window_size_4_lambda = peak_length
                lambda_bg = self.lambda_bg/self.scan_window*window_size_4_lambda                
                if self.nolambda:
                    # skip local lambda
                    local_lambda = lambda_bg
                    tlambda_peak = float(peak_num_tags)/peak_length*window_size_4_lambda
                else:
                    left_peak = peak_start+self.shift_size # go to middle point of the first fragment
                    right_peak = peak_end-self.shift_size  # go to middle point of the last fragment
                    left_lregion = peak_summit-self.lregion/2
                    left_sregion = peak_summit-self.sregion/2
                    right_lregion = peak_summit+self.lregion/2
                    right_sregion = peak_summit+self.sregion/2
                    #(cnum_10k,cnum_5k,cnum_1k,cnum_peak) = (0,0,0,0)
                    #(tnum_10k,tnum_5k,tnum_1k,tnum_peak) = (0,0,0,0)
                    (cnum_sregion, cnum_lregion, cnum_peak, tnum_sregion, tnum_lregion, tnum_peak) = (0,0,0,0,0,0)
                    #smallest = min(left_peak,left_10k,left_5k,left_1k)
                    #largest = max(right_peak,right_10k,right_5k,right_1k)
        
                    while index_ctag < len_ctags:
                        if get_read_start(ctags[index_ctag]) < left_lregion:
                            # go to next control tag
                            index_ctag+=1
                        elif right_lregion < get_read_start(ctags[index_ctag])\
                            or index_ctag + 1 >= len_ctags:
                            # finalize and go to next peak region
                            flag_find_ctag_locally = False
                            index_ctag = prev_index_ctag 
                            break
                        else:
                            if not flag_find_ctag_locally:
                                flag_find_ctag_locally = True
                                prev_index_ctag = index_ctag
                            p = get_read_start(ctags[index_ctag])
                            prob = get_c_read_prob(ctags[index_ctag])
                            if left_peak <= p <= right_peak:
                            #if peak_start <= p <= peak_end:
                                cnum_peak += prob
                            if left_sregion <= p <= right_sregion:
                                cnum_sregion += prob
                                cnum_lregion += prob
                            else:
                                cnum_lregion += prob
                            index_ctag += 1 # go to next tag
                    
                    inds_in_peak = []
                    while index_ttag < len_ttags:
                        if get_read_start(ttags[index_ttag]) < left_lregion:
                            # go to next treatment tag
                            index_ttag+=1
                        elif right_lregion < get_read_start(ttags[index_ttag]) \
                            or index_ttag + 1 >= len_ttags:
                            # finalize and go to next peak region
                            flag_find_ttag_locally = False
                            index_ttag = prev_index_ttag 
                            break
                        else:
                            if not flag_find_ttag_locally:
                                flag_find_ttag_locally = True
                                prev_index_ttag = index_ttag
                            p = get_read_start(ttags[index_ttag])
                            prob = get_t_read_prob(ttags[index_ctag])
                            if left_peak <= p <= right_peak:
                            #if peak_start <= p <= peak_end:
                                inds_in_peak.append(index_ttag)
                                tnum_peak += prob
                            if left_sregion <= p <= right_sregion:
                                tnum_sregion += prob
                                tnum_lregion += prob
                            else:
                                tnum_lregion += prob
                            index_ttag += 1 # go to next tag
                    if chrom not in all_peak_inds:
                        all_peak_inds[chrom] = []
                    all_peak_inds[chrom].append(inds_in_peak)
        
                    clambda_peak = float(cnum_peak)/peak_length*treat2control_ratio*window_size_4_lambda
                    #clambda_10k = float(cnum_10k)/self.third_lambda_region*treat2control_ratio*window_size_4_lambda
                    clambda_lregion = float(cnum_lregion)/self.lregion*treat2control_ratio*window_size_4_lambda
                    clambda_sregion = float(cnum_sregion)/self.sregion*treat2control_ratio*window_size_4_lambda
                    tlambda_peak = float(tnum_peak)/peak_length*window_size_4_lambda
                    #tlambda_10k = float(tnum_10k)/self.third_lambda_region*window_size_4_lambda
                    tlambda_lregion = float(tnum_lregion)/self.lregion*window_size_4_lambda
                    tlambda_sregion = float(tnum_sregion)/self.sregion*window_size_4_lambda
        
                    if pass_sregion:
                        # for experiment w/o control, peak region lambda and sregion region lambda are ignored!
                        local_lambda = max(lambda_bg,tlambda_lregion)
                    else:
                        # for experiment w/ control
                        outfile = open('lambdas.txt', 'a')
                        outfile.write('\t'.join(map(str, [lambda_bg,clambda_peak,clambda_lregion,clambda_sregion])) + '\n')
                        local_lambda = max(lambda_bg,clambda_peak,clambda_lregion,clambda_sregion)
                    if chrom not in all_local_lambdas:
                        all_local_lambdas[chrom] = []
                    all_local_lambdas[chrom].append(local_lambda)
        return all_peak_inds, all_local_lambdas
