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

class EffModelGC:
    def __init__(self, frags, prior, nr_cycles, mean_eff, max_eff, init_params, log, report):
        self.log        =   log
        self.report     =   report
        self.frags      =   frags
        self.prior      =   prior
        self.nr_cycles  =   nr_cycles
        self.init_params=   init_params
        self.mean_eff   =   mean_eff
        self.max_eff    =   max_eff

        self.calculate_ppr()
        self.fit()
        self.plot_report()

    def calculate_ppr(self):
        ratio           = np.zeros(self.frags.gc_marg.shape,dtype=float)
        z               = self.frags.frag_mat
        total           = float(self.frags.total_frags)
        size_marg       = self.frags.size_marg
        gc_size_prior   = self.prior.prior

        for gc in xrange(z.shape[1]):
            joint_tmp   = 0.0
            prior_tmp   = 0.0

            for size in xrange(z.shape[0]):
                prior       = size_marg[size] * gc_size_prior[size][gc]
                joint       = z[size, gc]/total
                joint_tmp   += joint
                prior_tmp   += prior

            if prior > 0.0:
                ratio[gc] = joint_tmp/prior_tmp

        self.ppr                = self.calculate_efficiencies(ratio)

    def calculate_efficiencies(self, ratios_vector):
        ratio                   = ratios_vector.copy()
        n                       = float(self.nr_cycles)
        Delta                   = self.calculate_delta(ratio)
        # Set the ratios < 1 to give zero efficiencies to avoid log problems:
        ratio[ratio == 0.0 ]    = np.exp( n * np.log(1/Delta) )
        # Calculate efficiencies: 
        ratio                   = np.exp(1.0/n * np.log(ratio) ) * Delta - 1
        # Discard negative "efficiencies":
        ratio[ratio < 0.0]      = 0.0
        return ratio

    def calculate_delta(self, ratio):
        Delta = None
        if self.max_eff != None:
            n                       = float(self.nr_cycles)
            max_ratio               = np.median(np.sort(ratio)[(len(ratio)-9):len(ratio)] )
            Delta                   = (1.0 + self.max_eff)/np.exp( 1.0/n * np.log(max_ratio) )
        elif self.mean_eff != None:
            Delta                   = (1.0 + self.mean_eff)
        else:
            self.log.fatal("No efficiency reference point specified!")
        return Delta
    
    def calculate_cost(self, params):
        # Calculate predicted efficiencies:
        gc = np.arange(0, 101,dtype=float)/100.0
        if len(self.ppr) != len(gc):
            self.log.fatal("Data vector has wrong size!")
        eff_pred    = self.gc_eff_func(gc, params)
        
        # Cost is weighted by p * (1- p)
        cost        = np.sum( ( (eff_pred - self.ppr) ** 2) * self.ppr * (-self.ppr + 1.0))
        return cost

    def fit(self):
        cost_func           = lambda p: self.calculate_cost(p)
        self.fitted_params  = optimize.fmin_tnc(func=cost_func, approx_grad=True, bounds=[(0,float('inf')), (0,1)], x0=self.init_params, disp=0)[0]
        self.gc_effs        = self.gc_eff_func(np.arange(0, 101,dtype=float)/100.0, self.fitted_params)

    def plot_report(self):
        fig                 = plt.figure()
        # Plot ppr:
        plt.bar(np.arange(len(self.ppr)),self.ppr,width=0.1)
        # Plot fitted curve:
        plt.plot(np.arange(len(self.gc_effs)),self.gc_effs)
        plt.title("Fitted GC efficiency function")
        self.report.pages.savefig(fig)

    def gc_eff_func(self,gc,params):
        a       =   float(params[0])
        b       =   float(params[1])
        effs    =   b + (1 - b) * np.power((1 - np.power(gc,a)), a)
        return effs

    def suggest_params(self):
        print "\t-eg \"(%g, %g, %g)\"" % (self.fitted_params[0], self.fitted_params[1], np.max(self.gc_effs))

    def save_raw(self, rj):
        rj.add_obj("nr_cycles",self.nr_cycles)
        y   = self.ppr
        x   = np.arange(len(y))
        rj.add_kv("gc_eff", x, y)
