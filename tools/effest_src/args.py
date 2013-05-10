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

VERSION = "1.1"

def parse_arguments():
    parser = argparse.ArgumentParser(description='Estimate GC-dependent fragment amplification efficiencies and fragment size distribution from paired-end RNA-seq data mapped to transcriptome (version %s).' % VERSION)

    parser.add_argument('input_file', metavar='input file', type=str, nargs='?',default=None,
                   help='Aligned *paired end* reads in SAM format sorted by *name*.')
    
    parser.add_argument('-f', metavar='ref_fasta', type=str, default=None,
                   help='Reference fasta.', required=False)
    
    parser.add_argument('-i', metavar='iso_list', type=str, default=None,
                   help='List of single isoform transcripts.', required=False)
   
    c_default   = 11 
    parser.add_argument('-c', metavar='nr_cycles', type=int, default=c_default,
                   help='Number of PCR cycles (%s).' % c_default, required=False)
   
    m_default   = 0.87
    parser.add_argument('-m', metavar='mean_eff', type=float, default=m_default,
                   help='Assumed pool efficiency (%s).' % m_default, required=False)

    parser.add_argument('-M', metavar='max_eff', type=float, default=None,
                   help='Assumed maximum efficiency (None).',required=False)

    parser.add_argument('-d', metavar='dist_fam', type=str, default="auto",
                   help='Distribution to model fragment size distribution (n|sn|*auto*).', required=False)

    parser.add_argument('-g', metavar='out_fasta', type=str, default=None,
                   help='Output fasta.', required=False)

    parser.add_argument('-j', metavar='out_json', type=str, default="raw_params.json",
                   help='File to store estimated raw parameters (raw_params.json).', required=False)

    e_default   = 10000.0
    parser.add_argument('-e', metavar='expr_mul', type=float, default=e_default,
                   help='Expression level multiplier (%s).' % e_default,required=False)

    u_default   = False
    parser.add_argument('-u', default=u_default, action="store_true", 
                   help='Save the uncorrected expression levels as well (False).', required=False)

    t_default   = True
    parser.add_argument('-t', default=t_default, action="store_true", 
                   help='Trim off old expression values (%s).' % t_default, required=False)
   
    w_default   = 5
    parser.add_argument('-w', metavar='step_size', type=int, default=w_default,
                   help='Sliding window size / step size ratio (%s).' % w_default, required=False)

    parser.add_argument('-s', metavar='out_count_file', type=str, default="effest_counts.pk",
                   help='Pickle counts to the specified file (effest_counts.pk).')

    parser.add_argument('-k', metavar='in_count_file', type=str, default=None,
                   help='Load counts from specifies file.')

    parser.add_argument('-p', metavar='out_prior_file', type=str, default="effest_pr.pk",
                   help='Pickle fragment prior to the specified file (effest_pr.pk).')

    parser.add_argument('-o', metavar='in_prior_file', type=str, default=None,
                   help='Load fragment prior from the specified pickle file.')

    q_default   = 0 
    parser.add_argument('-q', metavar='min_qual', type=int, default=q_default,
                   help='Minimum mapping quality (0).',required=False)

    parser.add_argument('-r', metavar='report_file', type=str, default="effest_report.pdf",
                   help='Report PDF (effest_report.pdf).')

    parser.add_argument('-l', metavar='log_file', type=str, default=None,
                   help='Log file.')
   
    v_default   = False
    parser.add_argument('-v', action='store_true' ,default=v_default,
                   help='Toggle verbose mode (%s).' % v_default)

    args                       =  parser.parse_args()
    if args.input_file is None and args.k is None:
        print >>sys.stderr, "No input alignment file or fragment pickle specified!"
        sys.exit(1)

    if args.r is None and (args.k is None or args.o is None):
        print >>sys.stderr, "Reference fasta not specified!"
        sys.exit(1)

    if args.m   == 0.0:
        args.m = None
    
    if args.M   == 0.0:
        args.M = None

    if args.m != None and args.M != None:
        args.m = None

    return args

def log_arguments(args, l):
    if args.m != None:
       l.vlog("Assumed mean efficiency: %g" % args.m) 
    
    if args.M != None:
       l.vlog("Assumed maximum efficiency: %g" % args.M) 

    l.vlog(  "Number of PCR cycles: %d"      % args.c)
    l.vlog( "Sliding window step ratio: %d" % args.w)
    l.vlog( "Minimum mapping quality: %d"   % args.q)
    l.vlog( "Saving estimated raw parameters to file: %s"   % args.j)

