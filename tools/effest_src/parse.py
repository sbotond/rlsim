#
# Copyright (C) 2013 EMBL - European Bioinformatics Institute
#
# This program is free software: you can redistribute it
# and/or modify it under the terms of the GNU General
# Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License
# for more details.
#
# Neither the institution name nor the name rlsim
# can be used to endorse or promote products derived from
# this software without prior written permission. For
# written permission, please contact <sbotond@ebi.ac.uk>.

# Products derived from this software may not be called
# rlsim nor may rlsim appear in their
# names without prior written permission of the developers.
# You should have received a copy of the GNU General Public
# License along with this program. If not, see
# <http://www.gnu.org/licenses/>.

import      HTSeq           as      ht
from        collections     import  defaultdict
import      string
import      math
import      warnings

class Frag:
    """ Class representing a fragment """
    def __init__(self, chrom, start, end, seq, strand, compact=False):
        self.chrom  =   chrom
        self.start  =   start
        self.end    =   end
        self.size   =   end - start
        self.strand =   "+"
        self.ttable =   string.maketrans("ATGC","TACG")

        if not compact:
            self.seq    =   seq
            if strand == "-":
                self.revcomp()

        self.gc     =   float( self.seq.count("G") + self.seq.count("C") )/float(self.size)
        self.seq    =   None

    def revcomp(self):
        """ reverse complement sequence """
        self.seq = self.seq.translate(self.ttable)[::-1]
        if self.strand == '+':
            self.strand = '-'
        else:
            self.strand = '+'

    def __str__(self):
        return ">%s %d -- %d %s (%d, %g)\n%s\n" % (self.chrom, self.start, self.end, self.strand, self.size, self.gc, self.seq)


class ParseFrags:
    """ Parse fragments """
    def __init__(self, infile, si, ref_seq, min_qual, log, report):
        self.si         =   si
        self.infile     =   infile
        self.min_qual   =   min_qual
        self.log        =   log
        self.report     =   report
        self.ref_seq    =   ref_seq
        self.read_iter  =   self.open_iter()

    def guess_input_type(self):
        in_type =   "SAM"
        fh      =   open(self.infile,"r")
        for c in fh.read(500):
            if c == "\0":
                in_type = "BAM"
                break
        fh.close()
        return in_type

    def open_iter(self):
        if self.infile == "-":
            self.infile = ht.FileOrSequence(sys.stdin)
        return ht.SAM_Reader(self.infile)

    def next_pair(self):
        """ Get next read pair """
        for (first, second) in ht.pair_SAM_alignments(self.read_iter):
            yield (first, second)

    def iter_frags(self):
        """ Get next fragment """
        min_qual    = self.min_qual
        for pair in self.next_pair():
            # Missing mate:
            if (pair[0] == None) or (pair[1] == None):
                continue
            # Unpaired read:
            elif (not pair[0].paired_end ) or (not pair[1].paired_end):
                self.log.fatal("Unpaired read found in alignment!")
            # Not a proper pair:
            elif (not pair[0].proper_pair ) or (not pair[1].proper_pair):
                continue
            # One of the reads is not aligned:
            elif (not pair[0].aligned ) or (not pair[1].aligned):
                continue
            elif (not pair[0].proper_pair ) or (not pair[1].proper_pair):
                continue
            # Mismatching reference:
            elif pair[0].iv.chrom != pair[1].iv.chrom:
                continue
            # One of the reads has low mapping quality:
            elif pair[0].aQual < min_qual or pair[1].aQual < min_qual:
                continue

            # Pull out useful info:
            chrom           =   pair[0].iv.chrom
            start, end      =   None, None
            strand          =   "+"

            # First read maps downstream:
            if  pair[0].iv.start < pair[1].iv.start:
                start       =   pair[0].iv.start
                end         =   pair[1].iv.end
            else:
            # First read maps upstream:
                strand      =   "-"
                start       =   pair[1].iv.start
                end         =   pair[0].iv.end

            # Get fragment sequence:
            frag_seq        =   self.ref_seq[chrom][start:end]
            yield Frag(chrom, start, end, frag_seq, strand)

    def parse_frags(self):
        """ Tabulate fragment statistics """
        infile_name = self.infile
        if type(self.infile) != str:
            infile_name = "<stdin>"
        self.log.vlog("Parsing fragments from file: %s" % infile_name)
        gc_mul      =   100
        frag_tab    =   defaultdict(lambda: defaultdict(float))
        # Dictionary to count fragments per transcript:
        expr_levels = { }

        # Catch and supress HTSeq warnings: 
        original_filters = warnings.filters[:]
        warnings.simplefilter("ignore")
        # Iterate over fragments:
        try:
            for frag in self.iter_frags():
                gc                          = int( math.floor(frag.gc * gc_mul) )
                # Only register fragments coming from single isoforms:
                if self.si.single_isoform(frag.chrom):
                        frag_tab[frag.size][gc]     += 1.0
                # Register all fragments to calculate fragment counts:
                if not expr_levels.has_key(frag.chrom):
                    expr_levels[frag.chrom] = np.zeros(101,dtype=float)
                tmp = expr_levels[frag.chrom]
                tmp[gc] += 1.0
        finally:
            # Restore the list of warning filters.
            warnings.filters = original_filters 

        # Transform dictionary into numpy array:
        max_size    =   max(frag_tab.keys())
        min_size    =   min(frag_tab.keys())
        self.min_fs =   min_size
        self.max_fs =   max_size

        frag_mat    =   np.zeros((max_size+1, 101),dtype=float)

        for (size, h) in frag_tab.iteritems():
            for (gc, count) in h.iteritems():
                frag_mat[size][gc]  = count

        frag_counts         = FragCounts(frag_mat, min_size, max_size, expr_levels, self.log, self.report)
        return frag_counts


class FragCounts:
    def __init__(self, frag_mat=None, min_size=None, max_size=None, expr_levels=None, log=None, report=None):
        if frag_mat == None:
            return
        self.log            = log
        self.report         = report
        self.frag_mat       = frag_mat
        self.min_size       = min_size
        self.max_size       = max_size
        self.expr_levels    = expr_levels
        self.total_frags    = np.sum(self.frag_mat)
        # Report the total number of fragments:
        note = ""
        self.log.log("Total number of %sfragments: %d" % (note, self.all_frags()))
        self.log.log("Total number of %sfragments from single isoform genes: %d" % (note, self.total_frags))
        # Calculate size marginal:
        self.size_marg  =   np.sum(self.frag_mat, axis=1)/float(self.total_frags)
        self.gc_marg    =   np.sum(self.frag_mat, axis=0)/float(self.total_frags)

    def all_frags(self):
        return sum( [np.sum(x) for x in self.expr_levels.values() ] )

    def plot_report(self):
        """ Plot various fragment statistics """
        self.report.plot_array(self.gc_marg, title="Marginal GC content distribution", xlab="GC content", ylab="Prob.")
        self.report.plot_array(self.size_marg, title="Marginal fragment size distribution",xlab="Fragment size", ylab="Prob.")
        self.report.plot_contour(self.frag_mat, title="Joint fragment distribution", xlab="GC", ylab="size")
                                                                                                
    def save_raw(self, rj):
        rj.add_obj("nr_frags", self.all_frags())
        y   = np.sum(self.frag_mat, axis=1) 
        x   = np.arange(len(y))
        rj.add_kv("frag_dist", x, y)

