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

import      sys
import      time
from        Bio                             import  SeqIO
from        matplotlib                      import  pyplot          as  plt
from        collections                     import  defaultdict
from        matplotlib.backends.backend_pdf import  PdfPages
import      numpy                           as      np

class Report:
    """ Class for plotting reports """
    def __init__(self, pdf):
        self.pdf    = pdf
        self.pages  = PdfPages(pdf)

    def plot_hash(self, h, title="", xlab="", ylab=""):
        """ Visualise hash as a bar plot """
        fig = plt.figure()
        plt.bar(h.keys(), h.values(), width=0.1)
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(title)
        self.pages.savefig(fig)
        plt.clf()
        plt.close(fig)

    def plot_contour(self, z, title="", xlab="", ylab=""):
        """ Visualise matrix as a filled contour plot """
        fig = plt.figure()
        p   = plt.contourf(z)
        plt.colorbar(p, orientation='vertical')
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(title)
        self.pages.savefig(fig)
        plt.clf()
        plt.close(fig)

    def plot_hist(self, data, title,xlab,ylab,bins): 
        fig = plt.figure()
        plt.hist(data, bins=bins)
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(title)
        self.pages.savefig(fig)
        plt.clf()
        plt.close(fig)

    def plot_array(self, y, title="", xlab="", ylab=""):
        """ Visualise  array as a bar plot """
        fig = plt.figure()
        plt.bar(np.arange(len(y)),y,width=0.1)
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(title)
        self.pages.savefig(fig)
        plt.clf()
        plt.close(fig)

    def close(self):
        self.pages.close()

class Log:
    """ Logging utility class """
    def __init__(self, fname=None, level=0):
        self.level = level
        if fname is None:
            self.fname  = "<sys.stderr>"     
            self.file   = sys.stderr
        else:
            self.file   = open(fname, "w")
            self.fname  = fname

    def close(self):
        self.file.flush()
        self.file.close()

    def log(self, message):
        if self.level < 0:
            return
        self.file.write("[%s] %s\n" % (time.strftime("%y-%m-%d %H:%M:%S"), message) )

    def vlog(self, message):
        if self.level < 1:
            return
        self.file.write("[%s] %s\n" % (time.strftime("%y-%m-%d %H:%M:%S"), message) )
         

    def fatal(self, message):
        self.file.write("[%s] %s\n" % (time.strftime("%y-%m-%d %H:%M:%S"), message) )
        sys.exit(1)

class Fasta:
    """ Fasta parsing class """
    def __init__(self, infile):
        self.infile     = infile
        self.in_fh      = open(infile, "r")
        self.iter       = SeqIO.parse(self.in_fh,'fasta')

    def __iter__(self):
        """ Return iterator """
        return iter(self.iter)

    def slurp(self):
        """ Slurp sequences """
        records = { }
        for s in iter(self):
            records[s.name] = (str(s.seq)).upper()
        return records

class SingleIsoList:
    """ The list of single isoform genes. """
    def __init__(self, infile, log):
        self.log        = log
        self.have_list  = False
        if not infile is None:
            self.infile = infile
            fh          = open(infile, "r")
            self.trs    = fh.readlines()
            for i in xrange(len(self.trs)):
                self.trs[i] = self.trs[i].rstrip()
            self.have_list = True
            self.log.vlog("Number of single isoform transcripts specified: %d" % self.size() )
        

    def single_isoform(self, name):
        if not self.have_list:
            return True
        else:
            return (name in self.trs) 

    def size(self):
        return len(self.trs)
