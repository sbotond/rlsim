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

class FragUnpickle:
    def __init__(self, infile, log, report):
        self.log    = log
        self.report = report
        self.infile = infile
        self.pickle = None
        if infile != None:
            self.pickle = cPickle.Unpickler( open(infile, "r") )

    def load(self):
        if not self.pickle is None:
            self.log.vlog("Loading fragment counts from file: %s" % self.infile)
            tmp = self.pickle.load()
            f   = FragCounts()
            f.__dict__ = tmp
            f.log       = self.log
            f.report    = self.report
            return f
        return None

class FragPickle:
    def __init__(self, outfile, log):
        self.log     = log
        self.pickle  = None
        self.outfile = outfile
        if outfile != None:
            self.pickle = cPickle.Pickler( open(outfile, "w") )

    def save(self, frags):
        self.log.vlog("Saving fragment counts to file: %s" % self.outfile)
        f   = { 
            "frag_mat": frags.frag_mat,
            "min_size": frags.min_size,
            "max_size": frags.max_size,
            "expr_levels": frags.expr_levels,
            "total_frags": frags.total_frags,
            "size_marg": frags.size_marg,
            "gc_marg": frags.gc_marg
        }
        if self.pickle:
            self.pickle.dump(f)

class PriorUnpickle:
    def __init__(self, infile, log, report):
        self.log    = log
        self.report = report
        self.pickle = None
        self.infile = infile
        if infile != None:
            self.pickle = cPickle.Unpickler( open(infile, "r") )

    def load(self):
        if self.pickle != None:
            self.log.vlog("Loading fragment prior from file: %s" % self.infile)
            tmp = self.pickle.load()
            f   = FragCounts()
            f.__dict__ = tmp
            f.log       = self.log
            f.report    = self.report
            return f
        return None

class PriorPickle:
    def __init__(self, outfile, log):
        self.log        = log
        self.pickle     = None
        self.outfile    = outfile
        if outfile != None:
            self.pickle = cPickle.Pickler( open(outfile, "w") )

    def save(self, pr):
        self.log.vlog("Saving fragment prior to file: %s" % self.outfile)
        f   = { 
            "prior": pr.prior,
            "gc_marg": pr.gc_marg,
        }
        if self.pickle:
            self.pickle.dump(f)
