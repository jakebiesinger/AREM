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
    * Added align_by_EM, incorporated it into call_peaks_w_control
    * Added diagnostic plots for EM steps
    
"""

import os
from math import log10 as math_log10
from array import array
from itertools import count as itertools_count, izip as itertools_izip

from AREM.OutputWriter import zwig_write
from AREM.IO.FeatIO import PeakIO,WigTrackI,BinKeeperI
from AREM.Prob import poisson_cdf,poisson_cdf_inv
from AREM.Constants import *

class PeakDetect:
    """Class to do the peak calling.

    e.g:
    >>> from AREM.PeakDetect import PeakDetect
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

    def call_peaks (self):
        """Call peaks function.

        Scan the whole genome for peaks. RESULTS WILL BE SAVED IN
        self.final_peaks and self.final_negative_peaks.
        """
        if self.control:                # w/ control
            self.peaks = self._call_peaks_w_control ()
        else:                           # w/o control
            self.peaks = self._call_peaks_wo_control ()
        return None

    def diag_result (self):
        """Run the diagnosis process on sequencing saturation.
        
        """
        if not self.diag:
            return None
        if self.control:                # w/ control
            return self._diag_w_control()
        else:                           # w/o control
            return self._diag_wo_control()

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
                text += "\t%.2f\t%.2f" % (peak[5],peak[6])
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
                text += "\t%.2f\t%.2f" % (peak[5],peak[6])
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
                text+= "%s\t%d\t%d\tAREM_peak_%d\t%.2f\n" % (chrom,peak[0],peak[1],n,peak[6])
        return text

    def summitsToBED (self):
        text = ""
        chrs = self.peaks.keys()
        chrs.sort()
        n = 0
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                n += 1
                text+= "%s\t%d\t%d\tAREM_peak_%d\t%.2f\n" % (chrom,peak[3]-1,peak[3],n,peak[4])
        return text

    def _add_fdr (self, final, negative): 
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

    def _call_peaks_w_control (self):
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
        self._shift_trackI(self.treat)
        self.info("#3 merge +/- strand of treatment data")

        self.treat.merge_plus_minus_locations_naive ()

        self.debug("#3 after shift and merging, tags: %d" % (self.treat.total))
        if self.save_wig:
            self.info("#3 save the shifted and merged tag counts into wiggle file...")
            #build wigtrack
            #if self.save_wig:
            #    treatwig = self._build_wigtrackI(self.treat,space=self.opt.space)
            if self.opt.wigextend:
                zwig_write(self.treat,self.opt.wig_dir_tr,self.zwig_tr,self.opt.wigextend,log=self.info,space=self.opt.space,single=self.opt.singlewig)
            else:
                zwig_write(self.treat,self.opt.wig_dir_tr,self.zwig_tr,self.d,log=self.info,space=self.opt.space,single=self.opt.singlewig)
        self.info("#3 call treatment peak candidates")
        peak_candidates = self._call_peaks_from_trackI (self.treat)
        
        self.info("#3 shift control data")
        self.info("#3 merge +/- strand of control data")
        self._shift_trackI(self.control)
        self.control.merge_plus_minus_locations_naive ()

        self.debug("#3 after shift and merging, tags: %d" % (self.control.total))
        if self.save_wig:
            self.info("#3 save the shifted and merged tag counts into wiggle file...")
            #build wigtrack
            #if self.save_score:
            #    controlbkI = self._build_binKeeperI(self.control,space=self.opt.space)
            if self.opt.wigextend:
                zwig_write(self.control,self.opt.wig_dir_ctl,self.zwig_ctl,self.opt.wigextend,log=self.info,space=self.opt.space,single=self.opt.singlewig)
            else:
                zwig_write(self.control,self.opt.wig_dir_ctl,self.zwig_ctl,self.d,log=self.info,space=self.opt.space,single=self.opt.singlewig)

        self.info("#3 call negative peak candidates")
        negative_peak_candidates = self._call_peaks_from_trackI (self.control)
        
        if self.treat.total_multi > 0 and not self.opt.no_EM:
            self.info("#3.5 Perform EM on treatment multi reads")
            self._align_by_EM(self.treat, self.control, peak_candidates, self.ratio_treat2control,expName='treatment')
        
        self.info("#3 use control data to filter peak candidates...")
        self.final_peaks = self._filter_w_control(peak_candidates,self.treat,self.control, self.ratio_treat2control,fake_when_missing=True)
        
        self.info("#3 find negative peaks by swapping treat and control")
        if self.treat.total_multi > 0 and self.control.total_multi > 0 \
                and not self.opt.no_EM:
            self.info("#3.5 Reset treatment alignment probabilities")
            # temporarily undo EM by resetting to prior probs
            self.treat.prob_aligns, self.treat.prior_aligns = self.treat.prior_aligns, self.treat.prob_aligns
            self.info("#3. Perform EM on control multi reads")
            self._align_by_EM(self.control, self.treat,
                               negative_peak_candidates,
                               float(self.control.total) / self.treat.total,
                               expName='control')
        self.info("#3 use treat data to filter negative peak candidates...")
        self.final_negative_peaks = self._filter_w_control(negative_peak_candidates,self.control,self.treat, 1.0/self.ratio_treat2control,fake_when_missing=True)
        
        if self.treat.total_multi > 0 and self.control.total_multi > 0 \
                and not self.opt.no_EM:
            # put back EM'ed probabilities
            self.treat.prob_aligns, self.treat.prior_aligns = self.treat.prior_aligns, self.treat.prob_aligns
        
        return self._add_fdr (self.final_peaks, self.final_negative_peaks)

    def _call_peaks_wo_control (self):
        """To call peaks w/o control data.

        """
        self.lambda_bg = float(self.scan_window)*self.treat.total/self.gsize
        self.debug("#3 background lambda: %.2f " % (self.lambda_bg))
        self.min_tags = poisson_cdf_inv(1-pow(10,self.pvalue/-10),self.lambda_bg)+1
        self.debug("#3 min tags: %d" % (self.min_tags))

        self.info("#3 shift treatment data")
        self._shift_trackI(self.treat)
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
        peak_candidates = self._call_peaks_from_trackI (self.treat)
        self.info("#3 use self to calculate local lambda and  filter peak candidates...")
        self.final_peaks = self._filter_w_control(peak_candidates,self.treat,self.treat,1,pass_sregion=True)
        return self.final_peaks

    def _print_peak_info (self, peak_info):
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

    def _filter_w_control (self, peak_info, treatment, control, treat2control_ratio, pass_sregion=False, write2wig= False, fake_when_missing=False ):
        """Use control data to calculate several lambda values around
        1k, 5k and 10k region around peak summit. Choose the highest
        one as local lambda, then calculate p-value in poisson
        distribution.
        
        New: Try to refine peaks by looking for enriched subpeaks

        Return value type in this format:
        a dictionary
        key value : chromosome
        items : array of (peak_start,peak_end,peak_length,peak_summit,peak_height,peak_num_tags,peak_pvalue,peak_fold_enrichment)
        
        """
        final_peak_info = {}
        chrs = peak_info.keys()
        chrs.sort()
        total = 0
        t_prob_aligns = treatment.prob_aligns
        c_prob_aligns = control.prob_aligns
        scan_window = self.scan_window
        lmax = max
        lpoisson_cdf = poisson_cdf
        for chrom in chrs:
            self.debug("#3 Chromosome %s" % (chrom))
            n_chrom = 0
            final_peak_info[chrom] = []
            peak_list = peak_info[chrom]
            try:
                (ctags,tmp) = control.get_locations_by_chr(chrom)
                (ctags_ind, tmp) = control.get_indexes_by_chr(chrom)
            except:
                self.warn("Missing %s data, skip it..." % (chrom))
                if fake_when_missing:
                    ctags = [-1,]
                    ctags_ind = [0]
                    self.warn("Fake a tag at %s:%d" % (chrom,-1))
                    tmp=[]
                else:
                    continue
            try:
                (ttags,tmp) = treatment.get_locations_by_chr(chrom)
                (ttags_ind, tmp) = treatment.get_indexes_by_chr(chrom)
            except:
                self.warn("Missing %s data, skip it..." % (chrom))
                if fake_when_missing:
                    ttags = [-1,]
                    ttags_ind = [0]
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

                #window_size_4_lambda = min(self.first_lambda_region,lmax(peak_length,scan_window))
                window_size_4_lambda = lmax(peak_length,scan_window)
                lambda_bg = self.lambda_bg/scan_window*window_size_4_lambda                
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
                    #largest = lmax(right_peak,right_10k,right_5k,right_1k)

                    #print 'index_ctag: %s, ctags[j] %s' % (index_ctag, ctags[index_ctag])
                    while index_ctag < len_ctags:
                        if ctags[index_ctag] < left_lregion:
                            # go to next control tag
                            index_ctag+=1
                        elif right_lregion < ctags[index_ctag] \
                            or index_ctag + 1 >= len_ctags:
                            # finalize and go to next peak region
                            flag_find_ctag_locally = False
                            index_ctag = prev_index_ctag 
                            break
                        else:
                            if not flag_find_ctag_locally:
                                flag_find_ctag_locally = True
                                prev_index_ctag = index_ctag
                            p = ctags[index_ctag]
                            prob = c_prob_aligns[ctags_ind[index_ctag]]
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

                    inds_in_peak = array('i',[])
                    inds_in_peak_append = inds_in_peak.append
                    while index_ttag < len_ttags:
                        if ttags[index_ttag] < left_lregion:
                            # go to next treatment tag
                            index_ttag+=1
                        elif right_lregion < ttags[index_ttag] \
                            or index_ttag + 1 >= len_ttags:
                            # finalize and go to next peak region
                            flag_find_ttag_locally = False
                            index_ttag = prev_index_ttag
                            break
                        else:
                            if not flag_find_ttag_locally:
                                flag_find_ttag_locally = True
                                prev_index_ttag = index_ttag
                            p = ttags[index_ttag]
                            prob = t_prob_aligns[ttags_ind[index_ttag]]
                            if left_peak <= p <= right_peak:
                            #if peak_start <= p <= peak_end:  # Jake-- again, seems this would be more accurate...
                                tnum_peak += prob
                                tnum_peak_total += 1
                                inds_in_peak_append(index_ttag)
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
                        local_lambda = lmax(lambda_bg,tlambda_lregion)
                    else:
                        # for experiment w/ control
                        local_lambda = lmax(lambda_bg,clambda_peak,clambda_lregion,clambda_sregion)

                p_tmp = lpoisson_cdf(tlambda_peak,local_lambda,lower=False)
                if p_tmp <= 0:
                    peak_pvalue = 3100
                else:
                    peak_pvalue = math_log10(p_tmp) * -10
                
                #print 'counts:', cnum_peak, cnum_sregion, cnum_lregion
                #print 'lambdas:', clambda_peak, clambda_sregion, clambda_lregion, lambda_bg
                

                if not self.opt.no_greedy_caller:
                    # build up sub peaks from the candidate we are iterating over
                    # by greedily adding tags to the current subpeak.  To avoid
                    # local minima, we always look at least scanwindow bases away
                    # get reads within first scanwindow
                    cpr_tags = [ttags[inds_in_peak[0]]]
                    cpr_tags_append = cpr_tags.append
                    cpr_probs = [t_prob_aligns[ttags_ind[inds_in_peak[0]]]]
                    cpr_probs_append = cpr_probs.append
                    cpr_probs_extend = cpr_probs.extend
                    cpr_mass = cpr_probs[-1]
                    middle_probs = []
                    middle_probs_append = middle_probs.append
                    j = 1  # right-most index we are considering
                    
                    # start with a scan_window slice of reads
                    while j < len(inds_in_peak):
                        cur_dist = ttags[inds_in_peak[j]] - cpr_tags[0]
                        if cur_dist <= scan_window: # always look scan_window away
                            cpr_tags_append(ttags[inds_in_peak[j]])
                            cpr_probs_append(t_prob_aligns[ttags_ind[inds_in_peak[j]]])
                            cpr_mass += cpr_probs[-1]
                            j += 1
                        else:
                            break
                    cpr_width = cpr_tags[-1] - cpr_tags[0]
                    lambda_width = lmax(cpr_width, scan_window)
                    cpr_pval = lpoisson_cdf(cpr_mass,local_lambda * lambda_width / window_size_4_lambda,lower=False)
                    
                    # add additional reads to current peak if they improve the enrichment
                    # but always look at reads within scan_window bases away
                    while j < len(inds_in_peak):
                        test_posn = ttags[inds_in_peak[j]]
                        cur_dist = test_posn - cpr_tags[-1]
                        if cur_dist <= scan_window:
                            test_width = test_posn - cpr_tags[0]
                            lambda_width = lmax(test_width, scan_window)
                            test_mass = cpr_mass + sum(middle_probs) + t_prob_aligns[ttags_ind[inds_in_peak[j]]]
                            test_pval = lpoisson_cdf(test_mass,local_lambda * lambda_width / window_size_4_lambda,lower=False)
                            if test_pval < cpr_pval:
                                # enrichment improved-- add the tag
                                cpr_tags_append(ttags[inds_in_peak[j]])
                                cpr_probs_extend(middle_probs)
                                cpr_probs_append(t_prob_aligns[ttags_ind[inds_in_peak[j]]])
                                cpr_mass = test_mass
                                cpr_width = test_width
                                cpr_pval = test_pval
                                del middle_probs[:]
                                #print 'accepted'
                            else:
                                # ignore the tag for now, include it when
                                # considering adjacent reads
                                middle_probs_append(t_prob_aligns[ttags_ind[inds_in_peak[j]]])
                                #print 'denied'
                            j += 1

                        # the next read is too far away so call previous region a peak
                        else:
                            if cpr_pval <= 0:
                                cpr_pval = 3100
                            else:
                                cpr_pval = math_log10(cpr_pval) * -10
                            if cpr_pval > self.pvalue:
                                # passes p-value threshold set by user
                                p_start, p_end, p_length, p_summit, p_height = self._tags_call_peak(cpr_tags, cpr_probs)
                                cpr_enrich = p_height / local_lambda * window_size_4_lambda / self.d
                                final_peak_info[chrom].append((p_start, p_end, p_length, p_summit, p_height, cpr_mass, cpr_pval, cpr_enrich))
                                n_chrom += 1
                                total += 1

                            # reset cpr
                            del cpr_tags[:], cpr_probs[:], middle_probs[:]
                            cpr_tags_append(ttags[inds_in_peak[j]])
                            cpr_probs_append(t_prob_aligns[ttags_ind[inds_in_peak[j]]])
                            cpr_mass = cpr_probs[-1]
                            # grab a new scan_window slice
                            while j < len(inds_in_peak):
                                cur_dist = ttags[inds_in_peak[j]] - cpr_tags[0]
                                if cur_dist <= scan_window:
                                    cpr_tags_append(ttags[inds_in_peak[j]])
                                    cpr_probs_append(t_prob_aligns[ttags_ind[inds_in_peak[j]]])
                                    cpr_mass += cpr_probs[-1]
                                    j += 1
                                else:
                                    break
                            cpr_width = cpr_tags[-1] - cpr_tags[0]
                            if cpr_width > 0:
                                lambda_width = lmax(cpr_width, scan_window)
                                cpr_pval = lpoisson_cdf(cpr_mass,local_lambda * lambda_width / window_size_4_lambda,lower=False)
                            else:
                                cpr_pval = 1.
                        # clip last tag from the left if it improves enrichment
                        if len(cpr_tags) > 3:
                            test_width = cpr_tags[-1] - cpr_tags[1]
                            test_mass = cpr_mass - cpr_probs[0]
                            lambda_width = lmax(test_width, scan_window)
                            test_pval = lpoisson_cdf(test_mass,local_lambda * lambda_width / window_size_4_lambda,lower=False)
                            if test_pval < cpr_pval:
                                # enrichment improved-- add the tag
                                cpr_tags.pop(0)
                                cpr_probs.pop(0)
                                cpr_mass = test_mass
                                cpr_width = test_width
                                cpr_pval = test_pval
                        
                    if len(cpr_tags) > 1:  # call last reads
                        if cpr_pval <= 0:
                            cpr_pval = 3100
                        else:
                            cpr_pval = math_log10(cpr_pval) * -10
                        #print 'calling peak!', cpr_pval, '>', self.pvalue, ' ?'
                        if cpr_pval > self.pvalue:
                            p_start, p_end, p_length, p_summit, p_height = self._tags_call_peak(cpr_tags, cpr_probs)
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

    def _call_peaks_from_trackI (self, trackI):
        """Call peak candidates from trackI data. Using every tag as
        step and scan the self.scan_window region around the tag. If
        tag number is greater than self.min_tags, then the position is
        recorded.

        Return: data in this format. (peak_start,peak_end,peak_length,peak_summit,peak_height,peak_num_tags)
        """
        prob_aligns = trackI.prob_aligns
        peak_candidates = {}
        self.debug("#3 search peak condidates...")
        chrs = trackI.get_chr_names()
        total = 0
        min_tags = self.min_tags
        scan_window = self.scan_window
        for chrom in chrs:
            self.debug("#3 Chromosome %s" % (chrom))
            n_chrom = 0
            peak_candidates[chrom] = []
            peak_cand_append = peak_candidates[chrom].append
            (tags,tmp) = trackI.get_locations_by_chr(chrom)
            (tags_ind, tmp) = trackI.get_indexes_by_chr(chrom)
            len_t = len(tags)
            cpr_tags = []       # Candidate Peak Region tags
            cpr_tags_append = cpr_tags.append
            cpr_tags_pop = cpr_tags.pop
            cpr_tags.extend(tags[:min_tags-1])
            number_cpr_tags = min_tags-1
            p = min_tags-1 # Next Tag Index
            cpr_probs = [prob_aligns[tags_ind[i]] for i in xrange(min_tags-1)]
            cpr_probs_append = cpr_probs.append
            cpr_probs_pop = cpr_probs.pop
            #cpr_indices = range(p+1)
            while p < len_t:
                if number_cpr_tags >= min_tags:
                    if tags[p] - cpr_tags[-min_tags+1] <= scan_window:
                        # add next tag, if the new tag is less than self.scan_window away from previous no. self.min_tags tag
                        cpr_tags_append(tags[p])
                        cpr_probs_append(prob_aligns[tags_ind[p]])
                        #cpr_indices.append(p)
                        number_cpr_tags += 1
                        p+=1
                    else:
                        # candidate peak region is ready, call peak...
                        #(peak_start,peak_end,peak_length,peak_summit,peak_height) = self._tags_call_peak (cpr_tags, cpr_probs)
                        peak_values = self._tags_call_peak (cpr_tags, cpr_probs)
                        peak_prob = sum(cpr_probs)
                        peak_cand_append(peak_values +(peak_prob,))
                        del cpr_tags[:]  # reset, but keep same variable
                        del cpr_probs[:]
                        cpr_tags_append(tags[p])
                        cpr_probs_append(prob_aligns[tags_ind[p]])
                        #cpr_indices = [p]
                        number_cpr_tags = 1
                        total += 1
                        n_chrom += 1
                        p += 1
                else:
                    # add next tag, but if the first one in cpr_tags
                    # is more than self.scan_window away from this new
                    # tag, remove the first one from the list
                    if tags[p] - cpr_tags[0] >= scan_window:
                        cpr_tags_pop(0)
                        cpr_probs_pop(0)
                        number_cpr_tags -= 1
                        #cpr_indices.pop(0)
                    cpr_tags_append(tags[p])
                    cpr_probs_append(prob_aligns[tags_ind[p]])
                    #cpr_indices.append(p)
                    number_cpr_tags += 1
                    p+=1
            self.debug("#3 peak candidates: %d" % (n_chrom))
        self.debug("#3 Total number of candidates: %d" % (total))
        return self._remove_overlapping_peaks(peak_candidates)
                
    def _tags_call_peak (self, tags, probs ):
        """Project tags to a line. Then find the highest point.

        """
        d = self.d
        start = tags[0]-d/2
        end = tags[-1]+d/2       # +1 or not?
        region_length = end - start
        line= [0]*region_length
        for tag, prob in itertools_izip(tags, probs):
            t_start = tag
            tag_projected_start = t_start-start-d/2
            tag_projected_end = t_start-start+d/2
            for i in xrange(tag_projected_start,tag_projected_end):
                line[i] += prob
        tops = []
        top_height = 0
        for i in xrange(len(line)):
            if line[i] > top_height:
                top_height = line[i]
                tops = [i]
            elif line[i] == top_height:
                tops.append(i)
        peak_summit = tops[len(tops)/2]+start
        return (start,end,region_length,peak_summit,top_height)

    def _shift_trackI (self, trackI):
        """Shift trackI data to right (for plus strand) or left (for
        minus strand).

        """
        chrs = trackI.get_chr_names()
        shift_size = self.shift_size
        for chrom in chrs:
            plus_tags, minus_tags = trackI.get_locations_by_chr(chrom)
            # plus
            for i in xrange(len(plus_tags)):
                plus_tags[i] += shift_size
            # minus
            for i in xrange(len(minus_tags)):
                minus_tags[i] -= shift_size
    
    def _build_wigtrackI (self, trackI, space=10):
        """Shift trackI then build a wigTrackI object.
        
        """
        chrs = trackI.get_chr_names()
        wigtrack = WigTrackI()
        wigtrack.span = space
        d = self.d
        step = 10000000 + 2*d
        prob_aligns = trackI.prob_aligns
        
        for chrom in chrs:
            tags = trackI.get_locations_by_chr(chrom)[0]
            tags_ind = trackI.get_indexes_by_chr(chrom)[0]
            l = len(tags)
            window_counts = array(BYTE4,[0]*step)
            startp = -1*d
            endp   = startp+step
            index_tag = 0
            prob = prob_aligns[tags_ind[index_tag]]
            while index_tag<l:
                s = tags[index_tag]-d/2     # start of tag
                e = s+d                     # end of tag
            
                if e < endp:
                    # project tag to window_counts line
                    ps = s-startp # projection start
                    pe = ps+d     # projection end
                    for i in xrange(ps,pe):
                        window_counts[i] += prob
                    index_tag += 1
                    prob = prob_aligns[tags_ind[index_tag]]
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

    def _diag_w_control (self):
        # sample
        sample_peaks = {}
        for i in xrange(90,10,-10):
            self.info("#3 diag: sample %d%%" % i)
            sample_peaks[i]=self._diag_peakfinding_w_control_sample(float(i)/(i+10))
        return self._overlap (self.final_peaks, sample_peaks,top=90,bottom=10,step=-10)

    def _diag_peakfinding_w_control_sample (self, percent):
        self.treat.sample(percent) # because sampling is after
                                   # shifting, track.total is used
                                   # now.
        self.control.sample(percent)
        ratio_treat2control = float(self.treat.total)/self.control.total

        self.lambda_bg = float(self.scan_window)*self.treat.total/self.gsize # bug fixed...
        self.min_tags = poisson_cdf_inv(1-pow(10,self.pvalue/-10),self.lambda_bg)+1

        self.debug("#3 diag: after shift and merging, treat: %d, control: %d" % (self.treat.total,self.control.total))
        self.info("#3 diag: call peak candidates")
        peak_candidates = self._call_peaks_from_trackI (self.treat)

        self.info("#3 diag: call negative peak candidates")
        negative_peak_candidates = self._call_peaks_from_trackI (self.control)
        
        self.info("#3 diag: use control data to filter peak candidates...")
        final_peaks_percent = self._filter_w_control(peak_candidates,self.treat,self.control, ratio_treat2control)
        return final_peaks_percent
        
    def _diag_wo_control (self):
        # sample
        sample_peaks = {}
        for i in xrange(90,10,-10):
            self.info("#3 diag: sample %d%%" % i)
            sample_peaks[i]=self._diag_peakfinding_wo_control_sample(float(i)/(i+10))
        return self._overlap (self.final_peaks, sample_peaks,top=90,bottom=10,step=-10)

    def _diag_peakfinding_wo_control_sample (self, percent):

        self.lambda_bg = float(self.scan_window)*self.treat.total/self.gsize # bug fixed...
        self.min_tags = poisson_cdf_inv(1-pow(10,self.pvalue/-10),self.lambda_bg)+1

        self.treat.sample(percent)
        self.debug("#3 diag: after shift and merging, tags: %d" % (self.treat.total))
        self.info("#3 diag: call peak candidates")
        peak_candidates = self._call_peaks_from_trackI (self.treat)
        self.info("#3 diag: use self to calculate local lambda and  filter peak candidates...")
        final_peaks_percent = self._filter_w_control(peak_candidates,self.treat,self.treat,1,pass_sregion=True) # bug fixed...
        return final_peaks_percent

    def _overlap (self, gold_peaks, sample_peaks, top=90,bottom=10,step=-10):
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
            
            r = self._overlap_fe(gold_peaks, sample_peaks, fe_low, fe_up, top, bottom, step)
            if r:
                diag_result.append(r)
        return diag_result

    def _overlap_fe (self, gold_peaks, sample_peaks, fe_low, fe_up, top, bottom, step):
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


    def _remove_overlapping_peaks (self, peaks ):
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


    def _align_by_EM(self, treatment, control, init_regions,
                      treat2control_ratio, show_graphs=None, expName='treatment'):
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
        self.info('#3.5 Gathering background lambdas')
        all_peak_inds, all_peak_lambdas = self._get_all_peak_lambdas(init_regions,
            treatment, control, treat2control_ratio)
        # filter out non-multiread peaks and save only counts of unique reads + multi indices
        min_score = self.opt.min_score if self.opt.min_score is not None else 1e-3
        max_score = self.opt.max_score if self.opt.max_score is not None else 2e3
        final_regions = {}
        peak_posns = {}
        t_prob_aligns = treatment.prob_aligns
        t_prior_aligns = treatment.prior_aligns
        t_enrich_scores = treatment.enrich_scores
        t_group_starts = treatment.group_starts
        lmax = max
        lpoisson_cdf = poisson_cdf
        if show_graphs is None:
            show_graphs = self.opt.show_graphs  # use cmd line param unless overridden
        #if show_graphs:
        in_candidate = [] #[False] * len(treatment.prob_aligns)  # if peak is in cand region
        for chrom in all_peak_inds.keys():
            try:
                (ttags,tmp) = treatment.get_locations_by_chr(chrom)
                (ttags_ind,tmp) = treatment.get_indexes_by_chr(chrom)
            except:
                continue
            for i in xrange(len(all_peak_inds[chrom])):
                local_lambda = all_peak_lambdas[chrom][i]
                peak_inds = all_peak_inds[chrom][i]
                unique_count = 0
                multi_inds = []
                if chrom not in peak_posns:
                    peak_posns[chrom] = []
                peak_posns[chrom].append((ttags[peak_inds[0]],
                                          ttags[peak_inds[-1]]))
                for ind in peak_inds:
                    if ttags_ind[ind] != 0:  # a multi read
                        multi_inds.append(ttags_ind[ind])  # get the index into the prob array
                        #if show_graphs:
                        #in_candidate[ttags[ind][1]] = True
                    else:
                        unique_count += 1
                if len(multi_inds) == 0:
                    continue  # no multi-reads in this peak.  Skip it
                else:
                    if chrom not in final_regions:
                        final_regions[chrom] = []
                    final_regions[chrom].append((local_lambda, unique_count, multi_inds))
        print 'total_multi: ', treatment.total_multi
        #print 'total number of peaks in candidate regions:', sum(1 for inc in in_candidate if inc)
        
        # for each iteration
        #if False:
        if show_graphs:
            try:
                self._plot_EM_state(0, final_regions, peak_posns, in_candidate, expName=expName)
            except:
                print 'Plotting failed!'
        prev_entropy = None
        no_prior_prob_map = self.opt.no_prior_prob_map
        for iteration in itertools_count(1):  # until convergence
            cur_entropy = list(self._get_read_entropy(treatment))
            if prev_entropy is not None:  # check for convergence
                denom = sum(ent ** 2. for ent in cur_entropy)
                if denom == 0.0:
                    difference = 0.0
                else:
                    difference = sum([(cur_entropy[i] - prev_entropy[i])**2. 
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
                    multi_mass = sum(t_prob_aligns[i] for i in multi_inds)
                    pvalue = lpoisson_cdf(multi_mass + unique_mass, local_lambda,
                                         lower=False)
                    if pvalue <= 0:
                        score = max_score
                    else:
                        score = lmax(-math_log10(pvalue), min_score)
                        score = min(score, max_score)
                    for i in multi_inds:
                        t_enrich_scores[i] = score

            # normalize the alignment probabilities for all multireads
            for i in xrange(1, len(t_group_starts)):
                group_start = t_group_starts[i]
                if i < len(t_group_starts) - 1:
                    group_end = t_group_starts[i+1]
                else:
                    group_end = len(t_group_starts)
                group_range = range(group_start, group_end)
                if not no_prior_prob_map: # posterior = map prob * enrichment
                    enrich_vals = [0.] * (group_end - group_start)
                    enrich_total = 0.
                    for j in group_range:
                        enr_val = t_enrich_scores[j] * t_prior_aligns[j]
                        enrich_vals[j - group_start] = enr_val
                        enrich_total += enr_val
                    for j in group_range:
                        t_prob_aligns[j] = enrich_vals[j-group_start] / enrich_total
                else:
                    enrich_total = sum([t_enrich_scores[j] for j in group_range])
                    for j in group_range:
                        t_prob_aligns[j] = t_enrich_scores[j] / enrich_total
            if show_graphs:
                try:
                    self._plot_EM_state(iteration, final_regions, peak_posns, in_candidate,expName=expName)
                except:
                    pass
            prev_entropy = cur_entropy
            # rinse and repeat (until convergence)
    
    def _plot_EM_state(self, iteration, final_regions, peak_posns, in_candidate, expName, output_summary=False, state=[False]):
        '''Plot the current data state. These may include:
           Entropy Histogram, CDF of enrichment scores and alignment probabilities,
           ratio of FG to BG scores or probabilities
        '''
        if not state[0]:
            import matplotlib
            matplotlib.use('Agg')
            from matplotlib import pyplot
            state[0] = True
        if output_summary:
            with open(self.opt.name +'_'+expName+ '_EMpeaks_%s.txt' % iteration, 'w') as outfile:
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
        self._plot_entropy_hist(iteration, self.treat,expName)
        self._plot_enrichment_CDF(iteration, self.treat,expName)
        self._plot_probability_CDF(iteration, self.treat,expName)
        #self._plot_probability_mass_ratio(iteration, self.treat, self.control,expName)
        self._plot_max_probability(iteration, self.treat,expName)

    def _get_read_entropy(self, treatment, normed=True, minScore=None):
        'generator for entropy in read alignments'
        t_prob_aligns = treatment.prob_aligns
        t_enrich_score = treatment.enrich_scores
        for i in xrange(len(treatment.group_starts)):
            group_start = treatment.group_starts[i]
            if i < len(treatment.group_starts) - 1:
                group_end = treatment.group_starts[i+1]
            else:
                group_end = len(treatment.prob_aligns)
            group_range = range(group_start, group_end)
            if minScore is None:
                probs = [t_prob_aligns[j] for j in group_range]
            else:
                scores = [t_enrich_score[j] - minScore for j in group_range]
                scoreTotal = sum(scores)
                if scoreTotal < 0.:
                    scoreTotal = 0.
                # renormalize, or if no alignments had score > minScore, consider them uniform
                if scoreTotal > 1e-5:
                    probs = [score / scoreTotal for score in scores]
                else:
                    probs = [1./len(group_range)] * len(group_range)
            entropies = [p * math_log10(p) if p > 0. else 0. for p in probs]
            if normed:
                yield -sum(entropies) / math_log10(len(group_range))
            else:
                yield -sum(entropies)

    def _plot_entropy_hist(self, iteration, treatment, expName):
        from matplotlib import pyplot
        entropy = list(self._get_read_entropy(treatment))
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
        pyplot.savefig(self.opt.name +'_'+expName+ '_entropy_%s.png' % iteration)
        pyplot.close()

    def _plot_enrichment_CDF(self, iteration, treatment, expName):
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
        pyplot.savefig(self.opt.name +'_'+expName+ '_CDF_enrichment_%s.png' % iteration)
        pyplot.close()
    
    def _plot_probability_CDF(self, iteration, treatment, expName):
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
        pyplot.savefig(self.opt.name +'_'+expName+ '_CDF_prob_%s.png' % iteration)
        pyplot.close()
    
    def _plot_max_probability(self, iteration, treatment,expName):
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
        pyplot.savefig(self.opt.name +'_'+expName+ '_max_prob_%s.png' % iteration)
        pyplot.close()


    def _get_all_peak_lambdas(self, peak_info, treatment, control, treat2control_ratio):
        """
        from MACS: calculate the local (max) lambda for each peak.
        Also returns all tag indices within each peak
        
        """
        pass_sregion = False
        chroms = sorted(peak_info.keys())
        all_peak_inds = {}
        all_local_lambdas = {}
        t_prob_aligns = treatment.prob_aligns
        c_prob_aligns = control.prob_aligns

        for chrom in chroms:
            peak_list = peak_info[chrom]
            try:
                (ctags,tmp) = control.get_locations_by_chr(chrom)
                (ctags_ind, tmp) = control.get_indexes_by_chr(chrom)
            except:
                self.warn("Missing %s data, skip it..." % (chrom))
                ctags = [-1,]
                ctags_ind = [0]
                self.warn("Fake a tag at %s:%d" % (chrom,-1))
                tmp=[]
            try:
                (ttags,tmp) = treatment.get_locations_by_chr(chrom)
                (ttags_ind, tmp) = treatment.get_indexes_by_chr(chrom)
            except:
                self.warn("Missing %s data, skip it..." % (chrom))
                ttags = [-1,]
                ttags_ind = [0]
                self.warn("Fake a tag at %s:%d" % (chrom,-1))
                tmp=[]
                
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
                        if ctags[index_ctag] < left_lregion:
                            # go to next control tag
                            index_ctag+=1
                        elif right_lregion < ctags[index_ctag] \
                            or index_ctag + 1 >= len_ctags:
                            # finalize and go to next peak region
                            flag_find_ctag_locally = False
                            index_ctag = prev_index_ctag 
                            break
                        else:
                            if not flag_find_ctag_locally:
                                flag_find_ctag_locally = True
                                prev_index_ctag = index_ctag
                            p = ctags[index_ctag]
                            prob = c_prob_aligns[ctags_ind[index_ctag]]
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
                        if ttags[index_ttag] < left_lregion:
                            # go to next treatment tag
                            index_ttag+=1
                        elif right_lregion < ttags[index_ttag] \
                            or index_ttag + 1 >= len_ttags:
                            # finalize and go to next peak region
                            flag_find_ttag_locally = False
                            index_ttag = prev_index_ttag 
                            break
                        else:
                            if not flag_find_ttag_locally:
                                flag_find_ttag_locally = True
                                prev_index_ttag = index_ttag
                            p = ttags[index_ttag]
                            prob = t_prob_aligns[ttags_ind[index_ttag]]
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
                        #outfile = open('lambdas.txt', 'a')
                        #outfile.write('\t'.join(map(str, [lambda_bg,clambda_peak,clambda_lregion,clambda_sregion])) + '\n')
                        local_lambda = max(lambda_bg,clambda_peak,clambda_lregion,clambda_sregion)
                    if chrom not in all_local_lambdas:
                        all_local_lambdas[chrom] = []
                    all_local_lambdas[chrom].append(local_lambda)
        return all_peak_inds, all_local_lambdas
