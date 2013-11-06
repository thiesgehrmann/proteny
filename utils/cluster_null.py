from scipy import stats;
import random;
import bisect;
import numpy as np;

from ibidas.utils.util import debug_here;

###############################################################################

class cluster_null_clt:

  name   = 'cluster_null_clt';
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

  name   = 'cluster_null_score_shuff';
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
    nue_a = kwargs.pop('nue_a');
    nue_b = kwargs.pop('nue_b');

    nue = nue_a + nue_b;

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

class cluster_null_score_shuff_max:

  name   = 'cluster_null_score_shuff_max'
  scores = None;
  hits   = None;
  hitA   = None;
  hitB   = None;
  perm   = {};
  nperm  = {};

  #############################################################################

  def __init__(self, **kwargs):
    self.scores = kwargs.pop('scores');
    self.nperm  = kwargs.pop('nperm') if 'nperm' in kwargs else 20000;
    self.hits   = kwargs.pop('hits');
    self.hitA   = kwargs.pop('hitA');
    self.hitB   = kwargs.pop('hitB');
    perm        = {};
  #edef

  #############################################################################

  def permute(self, n, nue_a, nue_b):
    k = (n, nue_a, nue_b);
    
    if not(k in self.perm):
      S = [ 0 ] * self.nperm;
      for i in xrange(self.nperm):
        s = random.sample(self.hits, n + nue_a + nue_b);
        match_s = 2 * sum([x[12] for x in s[0:n]]);
        ue_a_s  = sum([self.hits[self.hitA[(x[2], x[3])]][12] for x in s[n:n+nue_a]])
        ue_b_s  = sum([self.hits[self.hitB[(x[5], x[6])]][12] for x in s[n+nue_a:]])
        S[i] = match_s - ue_a_s - ue_b_s;
      #efor
      self.perm[k] = sorted(S);
    #fi
    
    return self.perm[k];
  #edef

  #############################################################################

  def pvalue(self, **kwargs):
    score = kwargs.pop('score');
    n     = kwargs.pop('n');
    nue_a = kwargs.pop('nue_a');
    nue_b = kwargs.pop('nue_b');

    if nue_a + nue_b >= n:
      return 1.0;
    #fi
    
    p = self.permute(n, nue_a, nue_b);
    
    i = bisect.bisect_left(p, score);
    
    return 1.0 - float(i) / float(self.nperm);
  #edef

  #############################################################################

#eclass

###############################################################################
