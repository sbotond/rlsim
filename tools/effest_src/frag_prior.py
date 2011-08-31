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

# Vectorized function for transforming 
# strings into "GC arrays":
def is_gc(c):
    if (c == 'G') or (c == "C"):
        return 1
    return 0 
is_gc_v = np.vectorize(is_gc,otypes='')

class FragPrior:
    """ Class to calculate GC content priors for a range of fragment sizes """
    def __init__(self, frags, si, ref_seq, step , log, report):
        self.frags      =   frags
        self.si         =   si
        self.ref_seq    =   ref_seq
        self.log        =   log
        self.report     =   report
        self.step       =   step
        self.log.vlog("Calculating fragment prior")
        self.slide()
        self.plot_report()

    def slide(self):
        """ Calculate P(gc|size) """
        prior_counts   = np.zeros((self.frags.max_size+1,101),dtype=float)
        for (name, seq) in ref_seq.iteritems():
            expr_level = 0
            if self.frags.expr_levels.has_key(name):
               expr_level = np.sum(self.frags.expr_levels[name]) 
            # Skip if expression level is zero:
            if expr_level == 0.0:
                continue
            # Skip if not on the list of single isoform genes:
            if self.si.single_isoform(name):
                self.slide_transcript(seq, expr_level, prior_counts)

        # Calculate marginal GC distribution:
        gc_marg         = np.sum(prior_counts,axis=0)
        gc_marg         = gc_marg/np.sum(gc_marg)
        self.gc_marg    = gc_marg

        # Normalize counts:
        for i in xrange(prior_counts.shape[0]):
            s = np.sum(prior_counts[i,])
            if s != 0.0:
                prior_counts[i,] = prior_counts[i,]/s
        self.prior  = prior_counts
    
    def slide_transcript(self, seq, expr_level, prior_counts):
        npseq       = is_gc_v(np.array(list(seq),dtype='S1'))
        size_marg   = self.frags.size_marg
        
        for winsize in xrange(self.frags.min_size, self.frags.max_size+1):
            step = winsize/self.step + 1
            for i in xrange(0,len(npseq)-winsize, step):
                gc_content                          =   int( ( float(npseq[i:i+winsize].sum())/winsize ) * 100)
                prior_counts[winsize,gc_content]    +=  expr_level

    def plot_report(self):
        self.report.plot_contour(self.prior, title="GC content prior",xlab="GC content",ylab="Fragment size")

