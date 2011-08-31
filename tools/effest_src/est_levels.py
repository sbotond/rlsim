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

class EstLevels:
    """ Class for quick-and-dirty estimating expression levels """
    def __init__(self, ref, frags, nr_cycles, model, out_file, expr_mul, ref_file, trim, gc_correct, log, report): 
        self.log        = log
        self.ref        = ref
        self.nr_cycles  = nr_cycles
        self.counts     = frags.expr_levels
        self.model      = model
        self.trim       = trim
        self.expr_mul   = expr_mul
        self.gc_correct = gc_correct
        self.report     = report
        # Set outpuf file name if missing:
        if out_file is None:
            if ref_file is None:
                self.log.fatal("Cannot estimate expression levels as output fasta file name is not specified in any way!")
            else:
                out_file = os.path.basename(ref_file).split(".fas")[0] + "_expr.fas"
        if os.path.exists(out_file):
            self.log.fatal("Output fasta file %s already exists! Aborting!" % out_file) 
        self.out_file   = out_file

        self.estimate()
        self.plot_levels()

    def estimate(self):
        """ Estimate expression levels. """
        # Precalculate efficiency weights:
        w       = np.power( ( 1/(self.model.ppr + 1.0) ), self.nr_cycles ) 
        c       = self.counts
        levels  = { }
        lengths = [ ]
        mul     = self.expr_mul

        for name, seq in self.ref.iteritems():
            level = 0 
            length= float(len(seq))

            if c.has_key(name):
                if self.gc_correct:
                    level = np.sum((c[name]/length) * w) * mul
                else:
                    level = np.sum(c[name]/length) * mul

            levels[name] = np.ceil(level)

        self.levels = levels

    def save_seq(self):
        """ Save estimated expression levels as augmented Fasta file. """
        self.log.vlog("Saving estimated expression levels to file: %s" % self.out_file)
        fh = open(self.out_file,'w')
        for name, seq in self.ref.iteritems():
            print >>fh, ">%s$%d\n%s\n" % (self.trim_name(name), self.levels[name], seq)
        fh.close()

    def trim_name(self, name):
        """ Trim off old expression levels. """
        if not self.trim:
            return name
        return name.split('$')[0]

    def plot_levels(self):
        """ Plot the distribution of estimated expression levels. """
        d = np.array(self.levels.values())
        self.report.plot_hist(data=d, 
            title="Distribution of estimated expression levels",
            xlab="Expression levels",
            ylab="Counts",
            bins=400
        )

        d = np.array(self.levels.values())
        d = d[d>0]
        self.report.plot_hist(data=d, 
            title="Distribution of estimated expression levels (zero excluded)",
            xlab="Expression levels",
            ylab="Counts",
            bins=400
        )
