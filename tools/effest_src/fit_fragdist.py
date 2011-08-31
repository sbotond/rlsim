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

class FitNormal:
    """ Class for fitting truncated normal distribution of fragment sizes. """
    def __init__(self, frags, dist_family, log, report):
        self.frags          = frags
        self.dist_family    = dist_family
        self.log            = log
        self.report         = report
        self.low            = self.frags.min_size
        self.high           = self.frags.max_size

    def tnorm_pdf(self,x, mean, sd, low, high):
        """ Calculate truncated normal density. """
        # Float conversions necessary in order ot avoid 
        # integer division.
        
        sd      = float(sd)
        mean    = float(mean)
        low     = float(low)
        high    = float(high)
        a       = (low - mean)/sd
        b       = (high -mean)/sd
        return stats.truncnorm.pdf(x, a, b, loc=mean, scale=sd)

    def tnorm_log_pdf(self,x, mean, sd, low, high):
        """ Calculate truncated normal density. """
        # Float conversions necessary in order ot avoid 
        # integer division.
        
        sd      = float(sd)
        mean    = float(mean)
        low     = float(low)
        high    = float(high)
        a       = (low - mean)/sd
        b       = (high -mean)/sd
        return stats.truncnorm.logpdf(x, a, b, loc=mean, scale=sd)


    def init_params(self):
        """ Initialize parameters based on sample mean and sd. """
        fd  = self.frags.size_marg

        m   = 0.0 
        for i in xrange(self.low, self.high+1):
            m += fd[i] * i

        s   = 0.0
        for i in xrange(self.low, self.high+1):
            s += fd[i] * (i - m) ** 2
        s = np.sqrt(s)
        return np.array([m,s],dtype=float) 

    def log_lik(self, p, x, d):
        """ Calculate log likelihood """
        pdf = self.tnorm_log_pdf(x, p[0], p[1], self.low, self.high)
        i   = pdf > float('-inf')
        L   = np.sum(pdf[i] * d[i])
        return L
        
    def fit(self):
        """ Fit truncated normal distribution """
        self.log.vlog("Estimating fragment size distribution (truncated normal)") 
        f           = self.frags.size_marg
        init_params = self.init_params()
        
        x       = np.arange(self.low, self.high+1)
        d       = np.sum(self.frags.frag_mat, axis=1)[self.low:self.high+1]

        cost_func           = lambda p: -self.log_lik(p,x,d)
        fitted_params   = optimize.fmin(func=cost_func, x0=init_params, maxiter=10**6, maxfun=10**6, disp=False)
        #fitted_params   = optimize.fmin(func=cost_func, x0=init_params, maxiter=10**6, maxfun=10**6, disp=True)
        #fitted_params  = optimize.fmin_tnc(func=cost_func, approx_grad=True, bounds=[(0,float('inf')), (0,float('inf'))], x0=init_params, disp=0)[0]
        self.L          = self.log_lik(fitted_params, x, d)

        self.mean   = fitted_params[0]
        self.sd     = fitted_params[1]

        self.log.vlog("Estimated mean fragment length: %d" % self.mean)
        self.log.vlog("Estimated fragment length standard deviation: %d" % self.sd)

    def report_model(self):
        """ Suggest parameters based on the fitted model. """
        print "Suggested rlsim parameters:"
        print "\t-d \"1.0:n:(%d,%d,%d,%d)\"" % (self.mean, self.sd, self.low, self.high)

    def plot_model(self):
        """ Plot fitted truncated normal. """ 
        fig                 = plt.figure()
        f                   = self.frags.size_marg/np.sum(self.frags.size_marg)
        # Plot empirical distribution:
        plt.bar(np.arange(len(f)),f,width=0.1)
        # Plot fitted curve:
        y = self.tnorm_pdf(x=np.arange(len(f)+1), low=self.low, high=self.high, mean=self.mean, sd=self.sd) 
        y[0:self.low] = 0.0
        plt.plot(np.arange(len(y)),y)
        plt.xlim((self.low-100, self.high+100))
        plt.title("Fitted fragment size distribution (truncated normal)") 
        self.report.pages.savefig(fig)

class FitSkewNormal:
    """ Fit skew normal distribution on the fragment size distribution. """
    def __init__(self, frags, dist_family, log, report):
        self.frags          = frags
        self.dist_family    = dist_family
        self.log            = log
        self.report         = report
        self.low            = self.frags.min_size
        self.high           = self.frags.max_size

    def snorm_pdf(self, x, loc, scale, shape):
        """ Calculate skew normal density. """
        x   = (x - loc)/float(scale)
        return 2.0 / scale * stats.norm.pdf(x) * stats.norm.cdf(shape*x)
    
    def init_params(self):
        """ Initialize parameters based on sample mean and sd. """
        fd  = self.frags.size_marg

        m   = 0.0 
        for i in xrange(self.low, self.high+1):
            m += fd[i] * i

        s   = 0.0
        for i in xrange(self.low, self.high+1):
            s += fd[i] * (i - m) ** 2
        s = np.sqrt(s)
        return np.array([m,s,1.0],dtype=float) 

    def log_lik(self, p, x, d):
        """ Calculate log likelihood """
        pdf = self.snorm_pdf(x, p[0], p[1], p[2])
        i   = pdf > 0.0
        L   = np.sum(np.log(pdf[i]) * d[i])
        return L
        
    def fit(self):
        """ Fit skew normal distribution. """
        self.log.vlog("Estimating fragment size distribution (skew normal)") 
        f           = self.frags.size_marg
        init_params = self.init_params()
        
        x       = np.arange(self.low, self.high+1)
        d       = np.sum(self.frags.frag_mat, axis=1)[self.low:self.high+1]

        cost_func           = lambda p: -self.log_lik(p,x,d)
        fitted_params  = optimize.fmin(func=cost_func, x0=init_params, maxiter=10**6, maxfun=10**6, disp=False)
        #fitted_params  = optimize.fmin(func=cost_func, x0=init_params, maxiter=10**6, maxfun=10**6, disp=True)
        #fitted_params  = optimize.fmin_tnc(func=cost_func, approx_grad=True, bounds=[ (0,float('inf')), (0,float('inf')), (float('-inf'),float('inf'))], x0=init_params, disp=0)[0]
        self.L          = self.log_lik(fitted_params, x, d)

        self.loc        = fitted_params[0]
        self.scale      = fitted_params[1]
        self.shape      = fitted_params[2]

        self.log.vlog("Estimated location parameter: %d" % self.loc)
        self.log.vlog("Estimated scale parameter: %d" % self.scale)
        self.log.vlog("Estimated shape parameter: %g" % self.shape)

    def report_model(self):
        """ Suggest parameters based on the fitted model. """
        print "Suggested rlsim parameters:"
        print "\t-d \"1.0:sn:(%d, %d, %g, %d, %d)\"" % (self.loc, self.scale, self.shape, self.low, self.high)

    def plot_model(self):
        """ Plot fitted skew normal distribution. """
        fig                 = plt.figure()
        f                   = self.frags.size_marg/np.sum(self.frags.size_marg)
        # Plot empirical distribution:
        plt.bar(np.arange(len(f)),f,width=0.1)
        # Plot fitted curve:
        #y = stats.truncnorm.pdf(np.arange(len(f)+1),self.low, self.high, loc=self.mean, scale=self.sd) 
        y = self.snorm_pdf(x=np.arange(self.high+1), loc=self.loc, scale=self.scale, shape=self.shape) 
        y = y/np.sum(y)
        plt.plot(np.arange(len(y)),y)
        plt.xlim((self.low-100, self.high+100))
        plt.title("Fitted fragment size distribution (skew normal)") 
        self.report.pages.savefig(fig)

def calc_aic(logL, nrp):
    return 2 * nrp - 2 * logL

def fit_fragdist(frags, dist_family, log, report):
    """ Fit a model on the fragment size distribution. """
    # Manually selected model:
    if dist_family == "n" or dist_family == "sn":
        if dist_family == "n":
            frag_model = FitNormal(frags, dist_family, log, report)
        elif dist_family == "sn":
            frag_model = FitSkewNormal(frags, dist_family, log, report)
        frag_model.fit()
        frag_model.report_model()
        frag_model.plot_model()

    # Model selection by AIC:
    elif dist_family == "auto":
        n    = FitNormal(frags, dist_family, log, report) 
        n.fit()
        n.plot_model()
        n_AIC   = calc_aic(n.L, 4) 
        L.vlog("Truncated normal AIC: %g" % n_AIC)

        sn   = FitSkewNormal(frags, dist_family, log, report) 
        sn.fit()
        sn.plot_model()
        sn_AIC  = calc_aic(sn.L, 3)
        L.vlog("Skew normal AIC: %g" % sn_AIC)
        
        L.vlog("Absolute AIC difference: %g" % abs(sn_AIC - n_AIC) )

        if sn_AIC < n_AIC:
            sn.report_model()
        else:
            n.report_model()

    # Illegal model type:
    else:
        log.fatal("Invalid distribution type specified for fragment size distribution: %s" % dist_family) 

