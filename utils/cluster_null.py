from scipy import stats;
import numpy as np;

###############################################################################

class cluster_null_clt:

  scores = None;
  mu     = None;
  s2     = None;

  #############################################################################

  def __init__(self, scores, Ea, Eb):
    self.scores = scores;
    self.estmu(Ea, Eb);
    self.ests2(Ea, Eb);
  #edef

  #############################################################################

  def estmu(self, Ea, Eb):
    n = Ea * Eb;
    self.mu = sum(self.scores) / n;
  #edef

  #############################################################################

  def ests2(self, Ea, Eb):
    n = Ea * Eb;
    nonzeros = [ (s - self.mu)**2 for s in self.scores ];
    zeros    = (self.mu**2) * (n - len(self.scores));
    self.s2 = (1.0/n) * (sum(nonzeros) + zeros);
  #edef

  #############################################################################

  def pvalue(self, Ea, Eb, score):

    n = Ea * Eb;

    z = (score - (n * self.mu)) / (n * self.s2);

    p = 1 - stats.norm.cdf(z);

    return p;
  #edef

  #############################################################################

#eclass

###############################################################################
