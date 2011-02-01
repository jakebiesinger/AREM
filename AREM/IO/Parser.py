# Time-stamp: <2011-01-20 18:21:42 Jake Biesinger>

"""Module for all AREM Parser classes for input.

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

@author: Biesinger, W Jacob B
@contact: jake.biesinger@gmail.com


Changes to this file since original release of MACS 1.4 (summer wishes):
  December/January 2011
    * Updated names (AREM, not MACS14)
    * Added MultiRead parser as a base class using alignment qualities
    * Updated SAM, BAM, and Bowtie parsers to use multi-read parser    
"""


# ------------------------------------
# python modules
# ------------------------------------
import logging
import struct
import gzip
from random import randrange as random_range
from math import log as mathlog
from operator import mul as op_multipy
from AREM.Constants import *
from AREM.IO.FeatIO import FWTrackII
# ------------------------------------
# constants
# ------------------------------------
__version__ = "Parser $Revision$"
__author__ = "Tao Liu <taoliu@jimmy.harvard.edu>"
__doc__ = "All Parser classes"

# ------------------------------------
# Misc functions
# ------------------------------------

def guess_parser ( fhd ):
    parser_dict = {"BED":BEDParser,
                   "ELAND":ELANDResultParser,
                   "ELANDMULTI":ELANDMultiParser,
                   "ELANDEXPORT":ELANDExportParser,
                   "SAM":SAMParser,
                   "BAM":BAMParser,
                   "BOWTIE":BowtieParser
                   }
    order_list = ("BAM",
                  "BED",
                  "ELAND",
                  "ELANDMULTI",
                  "ELANDEXPORT",
                  "SAM",
                  "BOWTIE",
                  )
    
    for f in order_list:
        p = parser_dict[f](fhd)
        s = p.sniff()
        if s:
            logging.info("Detected format is: %s" % (f) )
            return p
    raise Exception("Can't detect format!")

# ------------------------------------
# Classes
# ------------------------------------
class StrandFormatError(Exception):
    """Exception about strand format error.

    Example:
    raise StrandFormatError('Must be F or R','X')
    """
    def __init__ (self, string, strand):
        self.strand = strand
        self.string = string
        
    def __str__ (self):
        return repr( "Strand information can not be recognized in this line: \"%s\",\"%s\"" % (self.string,self.strand) )

class BaseQualityError(Exception):
    pass
    
class GenericParser:
    """Generic Parser class.

    Inherit this to write your own parser.
    """
    def __init__ (self, fhd):
        self.fhd = fhd
        return

    def tsize(self):
        return

    def build_fwtrack (self, opt):
        return 

    def _fw_parse_line (self, thisline ):
        return

    def sniff (self):
        try:
            t = self.tsize()
        except:
            self.fhd.seek(0)
            return False
        else:
            if t<=10 or t>=10000:
                self.fhd.seek(0)
                return False
            else:
                self.fhd.seek(0)
                return t

class MultiReadParser:
    """Parser capable of handling reads with more than one mapping.
    
    Inherit from this class to allow your parser to handle multi-reads.
    You must provide _fw_parse_line, which must return a tuple of 
    (chromosome,fpos,strand,tagname,qualstr,mismatches) instead of just 
    (chromosome,fpos,strand) as the other parsers do.
    
    NOTE: input data must be sorted by tagname to handle multireads correctly.
    """
    
    def build_fwtrack (self, opt):
        """Build FWTrackII from all lines, return a FWTrackII object.

        Handle multi-reads here by building a probability and enrichment index
        or select only one alignment from each multi-read. 
        Initial alignment probabilities are set from read/mismatch qualities
        or from a uniform distribution.
        """
        fwtrack = FWTrackII()
        i = 0
        m = 0
        recent_tags = []
        random_select_one_multi = opt.random_select_one_multi
        no_multi_reads = opt.no_multi_reads
        min_score = opt.min_score
        prior_prob_snp = opt.prior_prob_snp
        no_prior_prob_map = opt.no_prior_prob_map
        if opt.qual_scale == 'auto':
            opt.qual_scale = self._guess_qual_scale()
        if opt.qual_scale == 'sanger+33':
            qual_offset = 33
        elif opt.qual_scale == 'illumina+64':
            qual_offset = 64
        group_starts_append = fwtrack.group_starts.append
        fwtrack_add_loc = fwtrack.add_loc
        match_probs = {} # {(1,30):p(match|phred=30), (0,30):p(mismatch|phred=30)}

        for grouplines in self._group_by_name(self.fhd):
            fwtrack.total+=1  # in ratios, only count reads, not total alignments
            if len(grouplines) == 1:
                # uniquely mapping reads
                i+=1
                if i == 1000000:
                    m += 1
                    logging.info(" %d alignments read." % (m*1000000))
                    i=0
                chromosome, fpos, strand, qualstr, mismatches = grouplines[0]
                fwtrack_add_loc(chromosome,fpos,strand,0) # 0'th index => unique
            else:
                if no_multi_reads:  # throw away multi-reads
                    fwtrack.total -= 1
                    continue
                elif random_select_one_multi:  # choose one alignment at random
                    i+=1
                    if i == 1000000:
                        m += 1
                        logging.info(" %d alignments read." % (m*1000000))
                        i=0
                    randline = grouplines[random_range(len(grouplines))]
                    chromosome,fpos,strand,qualstr,mismatches = randline
                    fwtrack_add_loc(chromosome,fpos,strand,0)
                else:  # use all alignments probabilistically
                    group_starts_append(fwtrack.total_multi + 1)  # starts at 1 (0 reserved for unique reads)
                    if no_prior_prob_map:
                        # don't use map quality; just assume uniform priors
                        for (chromosome,fpos,strand, qualstr,mismatches) in grouplines:
                            i+=1
                            if i == 1000000:
                                m += 1
                                logging.info(" %d alignments read." % (m*1000000))
                                i=0
                            fwtrack.total_multi += 1
                            fwtrack_add_loc(chromosome,fpos,strand,
                                            fwtrack.total_multi)
                        normed_probs = [1./len(grouplines)] * len(grouplines)
                    else:
                        # TODO: might want to be working in log space-- if many mismatches, we'll lose precision
                        qualstr = grouplines[0][3]  # all quality strings are shared across the group
                        group_total_prob = 0.
                        group_probs = []
                        group_probs_append = group_probs.append
                        for (chromosome,fpos,strand, qualstr, mismatches) in grouplines:
                            i+=1
                            if i == 1000000:
                                m += 1
                                logging.info(" %d alignments read." % (m*1000000))
                                i=0
                            fwtrack.total_multi += 1
                            fwtrack_add_loc(chromosome,fpos,strand,
                                            fwtrack.total_multi)
                            mismatches = set(mismatches)
                            read_prob = 1.
                            # P(SNP) = prior probability a SNP occurs at any base
                            # P(SE) = probability there was a sequencing error (from PHRED)
                            # _P(Map|SNP,SE)__MATCH__SNP__SE_
                            #       0           0     0    0    # can't map here without explanation
                            #       1           0     0    1
                            #       1           0     1    0
                            #       1           0     1    1
                            #       1           1     0    0
                            #       1           1     0    1
                            #       0           1     1    0    # wouldn't map here if SNP, but sequencer read reference
                            #       1           1     1    1
                            # we are interested in P(Mapping | Match), which is equivalent to:
                            # \Sum_{SNP \in {0,1}, SE \in {0,1}} p(SNP) * p(SE) * p(Map|SE,SNP), or:
                            # p(Map|match = 0):
                            #     p(SE) + p(SNP) + p(SE)*p(SNP)
                            # p(Map|match = 1):
                            #    1 - (p(SE) + p(SE)*p(SNP))
                            for b in xrange(len(qualstr)):
                                tup = (b in mismatches,qualstr[b])
                                if tup in match_probs:
                                    prob = match_probs[tup]
                                elif tup[0]:  # mismatch
                                    p_seq_error = 10. ** ((qualstr[b]-qual_offset)/-10.)
                                    prob = p_seq_error + prior_prob_snp + p_seq_error * prior_prob_snp
                                    match_probs[tup] = prob
                                else:  # match
                                    p_seq_error = 10. ** ((qualstr[b]-qual_offset)/-10.)
                                    prob = 1. - (p_seq_error + p_seq_error * prior_prob_snp)
                                    match_probs[tup] = prob
                                read_prob *= prob
                            # quick & dirty check-- only looking at last base
                            assert qualstr[b] >= qual_offset # Specified quality scale yielded a negative phred score!  You probably have the wrong PHRED scale!
                            assert 0.<=read_prob<=1.  # error with map qualities
                            #raise BaseQualityError("Specified quality scale yielded a negative phred score!  You probably have the wrong PHRED scale!")
                            group_probs_append(read_prob)
                            group_total_prob += read_prob
                        normed_probs = [p / group_total_prob for p in group_probs]
                    fwtrack.prob_aligns.extend(normed_probs)
                    fwtrack.prior_aligns.extend(normed_probs)
                    fwtrack.enrich_scores.extend([min_score] * len(grouplines))
        return fwtrack
    
    def _guess_qual_scale(self):
        """Guess which quality scale the file is using by reading 
        the first 100 reads.
        """
        qual_min, qual_max = None, None
        n = 0
        m = 0
        while n<100 and m<1000:
            m += 1
            thisline = self.fhd.readline()
            (chromosome,fpos,strand, tagname, qualstr, mismatches
                ) = self._fw_parse_line(thisline)
            if not fpos or not chromosome:
                continue
            if qual_min is None or qual_max is None:
                qual_min = min(qualstr)
                qual_max = max(qualstr)
            else:
                qual_min = min(qual_min, min(qualstr))
                qual_max = max(qual_max, max(qualstr))
            n += 1
        self.fhd.seek(0)
        if qual_min < 33:
            raise BaseQualityError("Unrecognized scale for read quality: min %s, max %s" % (qual_min, qual_max))
        elif qual_min <= qual_max <= 73:
            # likely Sanger quals (ascii 33 to 73) since illumina quals would
            # have to have terrible read qualities to fall in this category
            return 'sanger+33'
        elif 64 <= qual_min <= qual_max <= 104:
            # likely Illumina quals (ascii 64 to 104)
            return 'illumina+64'
        else:
            raise BaseQualityError("Unrecognized scale for read quality: min %s, max %s" % (qual_min, qual_max))

    
    def _group_by_name(self, filelines):
        """Group multi-reads together by their tag name. 
        Yields a list of alignments for each read 
        
        NOTE: This only works if the filelines themselves are sorted by tagname        
        """
        cur_tagname = None
        aligns = []
        for thisline in filelines:
            (chromosome,fpos,strand,tagname,qualstr,mismatches
                ) = self._fw_parse_line(thisline)
            if not fpos or not chromosome:
                print 'skipping', thisline
                continue
            if cur_tagname is None:
                cur_tagname = tagname  # first line
            elif cur_tagname != tagname:
                yield aligns  # ID changed, so report all alignments for the previous read
                cur_tagname = tagname
                aligns = []
            aligns.append([chromosome, fpos, strand, qualstr, mismatches])
        if cur_tagname is not None:
            yield aligns  # might be one last alignment-- reported when ID changes
    

class BEDParser(GenericParser):
    """File Parser Class for tabular File.

    """
    def __init__ (self,fhd):
        self.fhd = fhd

    def tsize (self):
        s = 0
        n = 0
        m = 0
        while n<10 and m<1000:
            m += 1
            thisline = self.fhd.readline()
            (chromosome,fpos,strand) = self._fw_parse_line(thisline)
            if not fpos or not chromosome:
                continue
            thisline = thisline.rstrip()
            thisfields = thisline.split()
            s += int(thisfields[2])-int(thisfields[1])
            n += 1
        self.fhd.seek(0)
        return int(s/n)

    def build_fwtrack (self, opt):
        """Build FWTrackII from all lines, return a FWTrackII object.

        Note: All locations will be merged (exclude the same
        location) then sorted after the track is built.

        If both_strand is True, it will store strand information in
        FWTrackII object.

        if do_merge is False, it will not merge the same location after
        the track is built.
        """
        fwtrack = FWTrackII()
        i = 0
        m = 0
        for thisline in self.fhd:
            (chromosome,fpos,strand) = self._fw_parse_line(thisline)
            i+=1
            if i == 1000000:
                m += 1
                logging.info(" %d" % (m*1000000))
                i=0
            if not fpos or not chromosome:
                continue
            fwtrack.add_loc(chromosome,fpos,strand)
        return fwtrack
    
    def _fw_parse_line (self, thisline ):
        thisline = thisline.rstrip()
        if not thisline or thisline[:5]=="track" or thisline[:7]=="browser" or thisline[0]=="#": return ("comment line",None,None)

        thisfields = thisline.split('\t')
        chromname = thisfields[0]
        try:
            chromname = chromname[:chromname.rindex(".fa")]
        except ValueError:
            pass

        if len(thisfields) < 6 : # default pos strand if no strand
                                 # info can be found
            return (chromname,
                    int(thisfields[1]),
                    0)
        else:
            if thisfields[5] == "+":
                return (chromname,
                        int(thisfields[1]),
                        0)
            elif thisfields[5] == "-":
                return (chromname,
                        int(thisfields[2]),
                        1)
            else:
                raise StrandFormatError(thisline,thisfields[5])

class ELANDResultParser(GenericParser):
    """File Parser Class for tabular File.

    """
    def __init__ (self,fhd):
        """
        """
        self.fhd = fhd

    def tsize (self):
        s = 0
        n = 0
        m = 0
        while n<10 and m<1000:
            m += 1
            thisline = self.fhd.readline()
            (chromosome,fpos,strand) = self._fw_parse_line(thisline)
            if not fpos or not chromosome:
                continue
            thisline = thisline.rstrip()
            thisfields = thisline.split()
            s += len(thisfields[1])
            n += 1
        self.fhd.seek(0)
        return int(s/n)

    def build_fwtrack (self, opt):
        """Build FWTrackII from all lines, return a FWTrackII object.

        """
        fwtrack = FWTrackII()
        i = 0
        m = 0
        for thisline in self.fhd:
            (chromosome,fpos,strand) = self._fw_parse_line(thisline)
            i+=1
            if i == 1000000:
                m += 1
                logging.info(" %d" % (m*1000000))
                i=0
            if not fpos or not chromosome:
                continue
            fwtrack.add_loc(chromosome,fpos,strand)
        return fwtrack
    
    def _fw_parse_line (self, thisline ):
        #if thisline.startswith("#") or thisline.startswith("track") or thisline.startswith("browser"): return ("comment line",None,None) # comment line is skipped
        thisline = thisline.rstrip()
        if not thisline: return ("blank",None,None)

        thisfields = thisline.split()
        thistaglength = len(thisfields[1])

        if len(thisfields) <= 6:
            return ("blank",None,None)

        try:
            chromname = thisfields[6]
            chromname = chromname[:chromname.rindex(".fa")]
        except ValueError:
            pass

        if thisfields[2] == "U0" or thisfields[2]=="U1" or thisfields[2]=="U2":
            strand = thisfields[8]
            if strand == "F":
                return (chromname,
                        int(thisfields[7])-1,
                        0)
            elif strand == "R":
                return (chromname,
                        int(thisfields[7])+thistaglength-1,
                        1)
            else:
                raise StrandFormatError(thisline,strand)
        else:
            return (None,None,None)

class ELANDMultiParser(GenericParser):
    """File Parser Class for ELAND multi File.

    Note this parser can only work for s_N_eland_multi.txt format.

    Each line of the output file contains the following fields: 
    1. Sequence name 
    2. Sequence 
    3. Either NM, QC, RM (as described above) or the following: 
    4. x:y:z where x, y, and z are the number of exact, single-error, and 2-error matches 
    found 
    5. Blank, if no matches found or if too many matches found, or the following: 
    BAC_plus_vector.fa:163022R1,170128F2,E_coli.fa:3909847R1 
    This says there are two matches to BAC_plus_vector.fa: one in the reverse direction 
    starting at position 160322 with one error, one in the forward direction starting at 
    position 170128 with two errors. There is also a single-error match to E_coli.fa.
    """
    def __init__ (self,fhd):
        """
        """
        self.fhd = fhd

    def tsize (self, strict = False):
        s = 0
        n = 0
        m = 0
        while n<10 and m<1000:
            m += 1
            thisline = self.fhd.readline()
            (chromosome,fpos,strand) = self._fw_parse_line(thisline)
            if not fpos or not chromosome:
                continue
            thisline = thisline.rstrip()
            thisfields = thisline.split()
            s += len(thisfields[1])
            n += 1
        self.fhd.seek(0)
        return int(s/n)

    def build_fwtrack (self, opt):
        """Build FWTrackII from all lines, return a FWTrackII object.

        Note only the unique match for a tag is kept.
        """
        fwtrack = FWTrackII()
        i = 0
        m = 0
        for thisline in self.fhd:
            (chromosome,fpos,strand) = self._fw_parse_line(thisline)
            i+=1
            if i == 1000000:
                m += 1
                logging.info(" %d" % (m*1000000))
                i=0
            if not fpos or not chromosome:
                continue
            fwtrack.add_loc(chromosome,fpos,strand)
        return fwtrack
    
    def _fw_parse_line (self, thisline ):
        if not thisline: return (None,None,None)
        thisline = thisline.rstrip()
        if not thisline: return ("blank",None,None)

        #if thisline[0] == "#": return ("comment line",None,None) # comment line is skipped
        thisfields = thisline.split()
        thistagname = thisfields[0]        # name of tag
        thistaglength = len(thisfields[1]) # length of tag

        if len(thisfields) < 4:
            return (None,None,None)
        else:
            thistaghits = sum(map(int,thisfields[2].split(':')))
            if thistaghits > 1:
                # multiple hits
                return (None,None,None)
            else:
                (chromname,pos) = thisfields[3].split(':')

                try:
                    chromname = chromname[:chromname.rindex(".fa")]
                except ValueError:
                    pass
                
                strand  = pos[-2]
                if strand == "F":
                    return (chromname,
                            int(pos[:-2])-1,
                            0)
                elif strand == "R":
                    return (chromname,
                            int(pos[:-2])+thistaglength-1,
                            1)
                else:
                    raise StrandFormatError(thisline,strand)


class ELANDExportParser(GenericParser):
    """File Parser Class for ELAND Export File.

    """
    def __init__ (self,fhd):
        self.fhd = fhd

    def tsize (self):
        s = 0
        n = 0
        m = 0
        while n<10 and m<1000:
            m += 1
            thisline = self.fhd.readline()
            (chromosome,fpos,strand) = self._fw_parse_line(thisline)
            if not fpos or not chromosome:
                continue
            thisline = thisline.rstrip()
            thisfields = thisline.split("\t")
            s += len(thisfields[8])
            n += 1
        self.fhd.seek(0)
        return int(s/n)

    def build_fwtrack (self, opt):
        """Build FWTrackII from all lines, return a FWTrackII object.

        Note only the unique match for a tag is kept.
        """
        fwtrack = FWTrackII()
        i = 0
        m = 0
        for thisline in self.fhd:
            (chromosome,fpos,strand) = self._fw_parse_line(thisline)
            i+=1
            if i == 1000000:
                m += 1
                logging.info(" %d" % (m*1000000))
                i=0
            if not fpos or not chromosome:
                continue
            fwtrack.add_loc(chromosome,fpos,strand)
        return fwtrack
    
    def _fw_parse_line (self, thisline ):
        #if thisline.startswith("#") : return ("comment line",None,None) # comment line is skipped
        thisline = thisline.rstrip()
        if not thisline: return ("blank",None,None)
    
        thisfields = thisline.split("\t")

        if len(thisfields) > 12 and thisfields[12]:
            thisname = ":".join(thisfields[0:6])
            thistaglength = len(thisfields[8])
            strand = thisfields[13]
            if strand == "F":
                return (thisfields[10],int(thisfields[12])-1,0)
            elif strand == "R":
                return (thisfields[10],int(thisfields[12])+thistaglength-1,1)
            else:
                raise StrandFormatError(thisline,strand)
        else:
            return (None,None,None)

class PairEndELANDMultiParser(GenericParser):
    """File Parser Class for two ELAND multi Files for Pair-End sequencing.

    Note this parser can only work for s_N_eland_multi.txt format.

    Each line of the output file contains the following fields: 
    1. Sequence name 
    2. Sequence 
    3. Either NM, QC, RM (as described above) or the following: 
    4. x:y:z where x, y, and z are the number of exact, single-error, and 2-error matches 
    found 
    5. Blank, if no matches found or if too many matches found, or the following: 
    BAC_plus_vector.fa:163022R1,170128F2,E_coli.fa:3909847R1 
    This says there are two matches to BAC_plus_vector.fa: one in the reverse direction 
    starting at position 160322 with one error, one in the forward direction starting at 
    position 170128 with two errors. There is also a single-error match to E_coli.fa.
    """
    def __init__ (self,lfhd,rfhd):
        self.lfhd = lfhd
        self.rfhd = rfhd        

    def tsize (self):
        s = 0
        n = 0
        m = 0
        while n<10 and m<1000:
            m += 1
            thisline = self.lfhd.readline()
            thisline = thisline.rstrip()
            if not thisline: continue
            thisfields = thisline.split("\t")
            s += len(thisfields[1])
            n += 1
        self.fhd.seek(0)
        return int(s/n)

    def build_fwtrack (self, opt, dist=200):
        """Build FWTrackII from all lines, return a FWTrackII object.

        lfhd: the filehandler for left tag file
        rfhd: the filehandler for right tag file
        dist: the best distance between two tags in a pair

        The score system for pairing two tags:

        score = abs(abs(rtag-ltag)-200)+error4lefttag+error4righttag

        the smaller score the better pairing. If the score for a
        pairing is bigger than 200, this pairing will be discarded.

        Note only the best pair is kept. If there are over two best
        pairings, this pair of left and right tags will be discarded.

        Note, the orders in left tag file and right tag file must
        match, i.e., the Nth left tag must has the same name as the
        Nth right tag.

        Note, remove comment lines beforehand.
        """
        fwtrack = FWTrackII()
        i = 0
        m = 0
        lnext = self.lfhd.next
        rnext = self.rfhd.next
        self.dist = dist
        try:
            while 1:
                lline = lnext()
                rline = rnext()
                (chromname,fpos,strand) = self._fw_parse_line(lline,rline)

                i+=1
                if i == 1000000:
                    m += 1
                    logging.info(" %d" % (m*1000000))
                    i=0
                if not fpos or not chromname:
                    continue

                try:
                    chromname = chromname[:chromname.rindex(".fa")]
                except ValueError:
                    pass
                
                fwtrack.add_loc(chromname,fpos,strand)

        except StopIteration:
            pass
        return fwtrack
    
    def _fw_parse_line (self, leftline, rightline ):
        # >HWI-EAS275_5:4:100:340:1199/1	GTGCTGGTGGAGAGGGCAAACCACATTGACATGCT	2:1:0	chrI.fa:15061365F0,15068562F0,chrIV.fa:4783988R1
        # >HWI-EAS275_5:4:100:340:1199/2	GGTGGTGTGTCCCCCTCTCCACCAGCACTGCGGCT	3:0:0	chrI.fa:15061451R0,15068648R0,15071742R0

        leftfields = leftline.split()
        lefttaglength = len(leftfields[1]) # length of tag
        rightfields = rightline.split()
        righttaglength = len(rightfields[1]) # length of tag

        if len(rightfields) < 4 or len(leftfields) < 4:
            # one of the tag cann't be mapped to genome
            return (None,None,None)
        else:
            lefthits = self._parse_line_to_dict(leftfields[3])
            righthits = self._parse_line_to_dict(rightfields[3])            
            parings = []

            for seqname in lefthits.keys():
                if not righthits.has_key(seqname):
                    continue
                else:
                    leftpses = lefthits[seqname] # pse=position+strand+error
                    rightpses = righthits[seqname]
                    for (lp,ls,le) in leftpses:
                        for (rp,rs,re) in rightpses:
                            # try to pair them
                            if ls == 'F':
                                if rs == 'R':
                                    score = abs(abs(rp-lp)-self.dist)+le+re
                                    if score < 200:
                                        #parings.append((score,seqname,int((lp+rp)/2),0) )
                                        parings.append((score,seqname,lp,0))
                                else:
                                    # strands don't match
                                    continue
                            else:
                                if rs == 'F':
                                    score = abs(abs(rp-lp)-self.dist)+le+re
                                    if score < 200:
                                        #parings.append((score,seqname,int((lp+rp)/2),1) )
                                        parings.append((score,seqname,lp,1))
                                else:
                                    # strands don't match
                                    continue
            if not parings:
                return (None,None,None)
            parings.sort()
            if len(parings)>1 and parings[0][0] == parings[1][0]:
                # >2 best paring, reject!
                return (None,None,None)
            else:
                return parings[0][1:]                                
                    
    def _parse_line_to_dict ( self, linestr ):
        items = linestr.split(',')
        hits = {}
        for item in items:
            if item.find(':') != -1:
                # a seqname section
                (n,pse) = item.split(":") # pse=position+strand+error
                try:
                    n = n[:n.rindex(".fa")]
                except ValueError:
                    pass

                hits[n]=[]
                try:
                    sindex = pse.rindex('F')
                except ValueError:
                    sindex = pse.rindex('R')
                p = int(pse[:sindex])
                s = pse[sindex]
                e = int(pse[sindex+1:])
                hits[n].append((p,s,e))
            else:
                # only pse section
                try:
                    sindex = pse.rindex('F')
                except ValueError:
                    sindex = pse.rindex('R')
                p = int(pse[:sindex])
                s = pse[sindex]
                e = int(pse[sindex+1:])
                hits[n].append((p,s,e))
        return hits

### Contributed by Davide, modified by Tao
class SAMParser(GenericParser):
    """File Parser Class for SAM File.

    Each line of the output file contains at least: 
    1. Sequence name 
    2. Bitwise flag
    3. Reference name
    4. 1-based leftmost position of clipped alignment
    5. Mapping quality
    6. CIGAR string
    7. Mate Reference Name
    8. 1-based leftmost Mate Position
    9. Inferred insert size
    10. Query sequence on the same strand as the reference
    11. Query quality
    
    The bitwise flag is made like this:
    dec	meaning
    ---	-------
    1	paired read
    2	proper pair
    4	query unmapped
    8	mate unmapped
    16	strand of the query (1 -> reverse)
    32	strand of the mate
    64	first read in pair
    128	second read in pair
    256	alignment is not primary
    512	does not pass quality check
    1024	PCR or optical duplicate
    """
    def __init__ (self,fhd):
        """
        """
        self.fhd = fhd

    def tsize (self):
        s = 0
        n = 0
        m = 0
        while n<10 and m<1000:
            m += 1
            thisline = self.fhd.readline()
            (chromosome,fpos,strand) = self._fw_parse_line(thisline)
            if not fpos or not chromosome:
                continue
            thisline = thisline.rstrip()
            thisfields = thisline.split("\t")
            s += len(thisfields[9])
            n += 1
        self.fhd.seek(0)
        return int(s/n)

    def build_fwtrack (self, opt):
        """Build FWTrackII from all lines, return a FWTrackII object.

        Note only the unique match for a tag is kept.
        """
        fwtrack = FWTrackII()
        i = 0
        m = 0
        for thisline in self.fhd:
            (chromosome,fpos,strand) = self._fw_parse_line(thisline)
            i+=1
            if i == 1000000:
                m += 1
                logging.info(" %d" % (m*1000000))
                i=0
            if not fpos or not chromosome:
                continue
            fwtrack.add_loc(chromosome,fpos,strand)
        return fwtrack
    
    def _fw_parse_line (self, thisline ):
        thisline = thisline.rstrip()
        if not thisline: return ("blank",None,None)
        if thisline[0]=="@": return ("comment line",None,None) # header line started with '@' is skipped
        thisfields = thisline.split()
        thistagname = thisfields[0]         # name of tag
        thisref = thisfields[2]
        bwflag = int(thisfields[1])
        if bwflag & 4 or bwflag & 512 or bwflag & 1024:
            return (None, None, None)       #unmapped sequence or bad sequence
        if bwflag & 1:
            # paired read. We should only keep sequence if the mate is mapped
            # and if this is the left mate, all is within  the flag! 
            if not bwflag & 2:
                return (None, None, None)   # not a proper pair
            if bwflag & 8:
                return (None, None, None)   # the mate is unmapped
            p1pos = int(thisfields[3]) - 1
            p2pos = int(thisfields[7]) - 1
            if p1pos > p2pos:
                # this pair is the farthest one, skip it
                return (None, None, None)
        # In case of paired-end we have now skipped all possible "bad" pairs
        # in case of proper pair we have skipped the rightmost one... if the leftmost pair comes
        # we can treat it as a single read, so just check the strand and calculate its
        # start position... hope I'm right!
        if bwflag & 16:
            thisstrand = 1
            thisstart = int(thisfields[3]) - 1 + len(thisfields[9])	#reverse strand should be shifted len(query) bp 
        else:
            thisstrand = 0
            thisstart = int(thisfields[3]) - 1	

        try:
            thisref = thisref[:thisref.rindex(".fa")]
        except ValueError:
            pass
        return (thisref, thisstart, thisstrand)

class BAMParser(GenericParser):
    """File Parser Class for BAM File.

    File is gzip-compatible and binary.
    Information available is the same that is in SAM format.
    
    The bitwise flag is made like this:
    dec	meaning
    ---	-------
    1	paired read
    2	proper pair
    4	query unmapped
    8	mate unmapped
    16	strand of the query (1 -> reverse)
    32	strand of the mate
    64	first read in pair
    128	second read in pair
    256	alignment is not primary
    512	does not pass quality check
    1024	PCR or optical duplicate
    """
    def __init__ (self,fhd):
        """
        """
        self.fhd = fhd

    def sniff(self):
        if self.fhd.read(3) == "BAM":
            return True
        else:
            return False
            
    def tsize(self):
        fseek = self.fhd.seek
        fread = self.fhd.read
        ftell = self.fhd.tell
        # move to pos 4, there starts something
        fseek(4)
        header_len =  struct.unpack('<i', fread(4))[0]
        fseek(header_len + ftell())
        # get the number of chromosome
        nc = struct.unpack('<i', fread(4))[0]
        #print nc
        for x in range(nc):
            # read each chromosome name
            nlength = struct.unpack('<i', fread(4))[0]
            # jump over chromosome size, we don't need it
            fread(nlength)
            fseek(ftell() + 4)
        s = 0
        n = 0
        while n<10:
            #print n
            entrylength = struct.unpack('<i', fread(4))[0]
            #print entrylength
            data = fread(entrylength)
            #print n
            #(chrid,fpos,strand) = self._fw_binary_parse(fread(entrylength))
            s += struct.unpack('<i', data[16:20])[0]
            #print n
            n += 1
            #print s,n
        fseek(0)
        logging.info("tag size: %d" % int(s/n))
        return int(s/n)

    def build_fwtrack (self, opt):
        """Build FWTrackII from all lines, return a FWTrackII object.

        Note only the unique match for a tag is kept.
        """
        fwtrack = FWTrackII()
        i = 0
        m = 0
        references = []
        fseek = self.fhd.seek
        fread = self.fhd.read
        ftell = self.fhd.tell
        # move to pos 4, there starts something
        fseek(4)
        header_len =  struct.unpack('<i', fread(4))[0]
        fseek(header_len + ftell())
        # get the number of chromosome
        nc = struct.unpack('<i', fread(4))[0]
        for x in range(nc):
            # read each chromosome name
            nlength = struct.unpack('<i', fread(4))[0]
            references.append(fread(nlength)[:-1])
            # jump over chromosome size, we don't need it
            fseek(ftell() + 4)
        
        while 1:
            try:
                entrylength = struct.unpack('<i', fread(4))[0]
            except struct.error:
                break
            (chrid,fpos,strand) = self._fw_binary_parse(fread(entrylength))                    
            i+=1
            if i == 1000000:
                m += 1
                logging.info(" %d" % (m*1000000))
                i=0
            if fpos >= 0:
                fwtrack.add_loc(references[chrid],fpos,strand)
        self.fhd.close()
        return fwtrack
    
    def _fw_binary_parse (self, data ):
        # we skip lot of the available information in data (i.e. tag name, quality etc etc)
        if not data: return (None,-1,None)

        thisref = struct.unpack('<i', data[0:4])[0]
        thisstart = struct.unpack('<i', data[4:8])[0]
        (cigar, bwflag) = struct.unpack('<hh', data[12:16])
        if bwflag & 4 or bwflag & 512 or bwflag & 1024:
            return (None, -1, None)       #unmapped sequence or bad sequence
        if bwflag & 1:
            # paired read. We should only keep sequence if the mate is mapped
            # and if this is the left mate, all is within  the flag! 
            if not bwflag & 2:
                return (None, -1, None)   # not a proper pair
            if bwflag & 8:
                return (None, -1, None)   # the mate is unmapped
            p1pos = thisstart
            p2pos = struct.unpack('<i', data[24:28])[0]
            if p1pos > p2pos:
                # this pair is the farthest one, skip it
                return (None, -1, None)
        # In case of paired-end we have now skipped all possible "bad" pairs
        # in case of proper pair we have skipped the rightmost one... if the leftmost pair comes
        # we can treat it as a single read, so just check the strand and calculate its
        # start position... hope I'm right!
        if bwflag & 16:
            thisstrand = 1
            thisstart = thisstart + struct.unpack('<i', data[16:20])[0]	#reverse strand should be shifted len(query) bp 
        else:
            thisstrand = 0

        return (thisref, thisstart, thisstrand)

### End ###

class BowtieParser(MultiReadParser, GenericParser):
    """File Parser Class for map files from Bowtie or MAQ's maqview
    program.

    """
    def __init__ (self,fhd):
        """
        """
        self.fhd = fhd

    def tsize (self):
        s = 0
        n = 0
        m = 0
        while n<10 and m<1000:
            m += 1
            thisline = self.fhd.readline()
            (chromosome,fpos,strand, tagname, qualstr, mismatches
                ) = self._fw_parse_line(thisline)
            if not fpos or not chromosome:
                continue
            thisfields = thisline.split('\t')
            s += len(thisfields[4])
            n += 1
        self.fhd.seek(0)
        return int(s/n)    
    
    def _fw_parse_line (self, thisline ):
        """
        The following definition comes from bowtie website:
        
        The bowtie aligner outputs each alignment on a separate
        line. Each line is a collection of 8 fields separated by tabs;
        from left to right, the fields are:

        1. Name of read that aligned

        2. Orientation of read in the alignment, - for reverse
        complement, + otherwise

        3. Name of reference sequence where alignment occurs, or
        ordinal ID if no name was provided

        4. 0-based offset into the forward reference strand where
        leftmost character of the alignment occurs

        5. Read sequence (reverse-complemented if orientation is -)

        6. ASCII-encoded read qualities (reversed if orientation is
        -). The encoded quality values are on the Phred scale and the
        encoding is ASCII-offset by 33 (ASCII char !).

        7. Number of other instances where the same read aligns
        against the same reference characters as were aligned against
        in this alignment. This is not the number of other places the
        read aligns with the same number of mismatches. The number in
        this column is generally not a good proxy for that number
        (e.g., the number in this column may be '0' while the number
        of other alignments with the same number of mismatches might
        be large). This column was previously described as"Reserved".

        8. Comma-separated list of mismatch descriptors. If there are
        no mismatches in the alignment, this field is empty. A single
        descriptor has the format offset:reference-base>read-base. The
        offset is expressed as a 0-based offset from the high-quality
        (5') end of the read.

        """
        thisline = thisline.rstrip('\n')
        if not thisline: return ("blank",None,None)
        if thisline[0]=="#": return ("comment line",None,None) # comment line is skipped
        thisfields = thisline.split('\t')
        
        tagname = thisfields[0]
        chromname = thisfields[2]
        qualstr = thisfields[5]
        qualstr = struct.unpack('%sb'%len(qualstr), qualstr)  # pretend the qualstr is a short array
        # mismatches of form: 2:T>G,38:C>A,42:C>A
        mismatches = thisfields[7]
        if len(mismatches) > 0:
            mismatches = [int(sub[:sub.index(':')]) for sub in 
                          mismatches.split(',')]
        else:
            mismatches = []
        try:
            chromname = chromname[:chromname.rindex(".fa")]
        except ValueError:
            pass

            if thisfields[1] == "+":
                return (chromname,
                        int(thisfields[3]),
                        0, tagname, qualstr, mismatches)
            elif thisfields[1] == "-":
                return (chromname,
                        int(thisfields[3])+len(thisfields[4]),
                        1, tagname, qualstr, mismatches)
            else:
                raise StrandFormatError(thisline,thisfields[1])

