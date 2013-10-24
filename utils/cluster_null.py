from scipy import stats;
import random;
import bisect;
import numpy as np;

from ibidas.utils.util import debug_here;

###############################################################################

class cluster_null_clt:

  scores = None;
  mu     = None;
  s2     = None;

  #############################################################################

  def __init__(self, **kwargs):
    self.scores = kwargs.pop('scores')
    Ea          = kwargs.pop('Ea');
    Eb          = kwargs.pop('Eb');
    self.estmu(Ea, Eb);
    self.ests2(Ea, Eb);
  #edef

  #############################################################################

  def estmu(self, Ea, Eb):
    n  = Ea * Eb;

    self.mu = sum(self.scores) / n;
  #edef

  #############################################################################

  def ests2(self, Ea, Eb):
    n  = Ea * Eb;

    nonzeros = [ (s - self.mu)**2 for s in self.scores ];
    zeros    = (self.mu**2) * (n - len(self.scores));
    self.s2 = (1.0/n) * (sum(nonzeros) + zeros);
  #edef

  #############################################################################

  def pvalue(self, **kwargs):
    score = kwargs.pop('score');
    Ea = kwargs.pop('Ea');
    Eb = kwargs.pop('Eb');
    n  = Ea * Eb;

    z = (score - (n * self.mu)) / (n * self.s2);

    p = 1 - stats.norm.cdf(z);

    return p;
  #edef

  #############################################################################

#eclass

###############################################################################

class cluster_null_score_shuff:

  scores = None;
  perm   = {};
  nperm  = {};

  #############################################################################

  def __init__(self, **kwargs):
    self.scores = kwargs.pop('scores');
    self.nperm  = kwargs.pop('nperm') if 'nperm' in kwargs else 20000;
    perm        = {};
  #edef

  #############################################################################

  def permute(self, n, nue):
    k = (n, nue);

    if not(k in self.perm):
      S = [ 0 ] * self.nperm;
      for i in xrange(self.nperm):
        s = random.sample(self.scores, n);
        S[i] = 2*sum(s[0:n]) - nue;
      #efor
      self.perm[k] = sorted(S);
    #fi

    return self.perm[k];
  #edef

  #############################################################################

  def pvalue(self, **kwargs):
    score = kwargs.pop('score');
    n     = kwargs.pop('n');
    nue   = kwargs.pop('nue');

    if nue >= n:
      return 1.0;
    #fi

    p = self.permute(n, nue);

    i = bisect.bisect_left(p, score);

    return 1.0 - float(i) / float(self.nperm);
  #edef

  #############################################################################

#eclass

###############################################################################
