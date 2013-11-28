from scipy import stats;
import random;
import bisect;
import numpy as np;

from ibidas.utils.util import debug_here;

###############################################################################

class cluster_null_clt:

  name   = 'ClusterNullCLT';
  H_mu   = None;
  H_s2   = None;
  U1_mu  = None;
  U1_s2  = None;
  U2_mu  = None;
  U2_s2  = None;

  #############################################################################

  def __init__(self, **kwargs):
    scores = kwargs.pop('scores');
    hits   = kwargs.pop('hits');
    hitA   = kwargs.pop('hitA');
    hitB   = kwargs.pop('hitB');    

    self.H_mu  = self.estmu(scores);
    self.H_s2  = self.ests2(scores);

    L          = [ hits[hitA[k]][12] for k in hitA.keys() ];
    self.U1_mu = self.estmu(L);
    self.U1_s2 = self.ests2(L);

    L          = [ hits[hitB[k]][12] for k in hitB.keys() ];
    self.U2_mu = self.estmu(L);
    self.U2_s2 = self.ests2(L);
  #edef

  #############################################################################

  def estmu(self, L):
    mu = np.mean(L);
    return mu;
  #edef

  #############################################################################

  def ests2(self, L):
    s2 = np.var(L);
    return s2;
  #edef

  #############################################################################

  def pvalue(self, **kwargs):
    # The score is a sum of CLT estimated normal distributions.

    score = kwargs.pop('score');
    n     = kwargs.pop('n');
    nue_a = kwargs.pop('nue_a');
    nue_b = kwargs.pop('nue_b');

    mu = (self.H_mu * n) - (self.U1_mu * nue_a) - (self.U2_mu * nue_b);
    s2 = (self.H_s2 * n) + (self.U1_s2 * nue_a) + (self.U2_s2 * nue_b);

    z = ( score - mu ) / s2;
    p = 1 - stats.norm.cdf(z);

    return p;
  #edef

  #############################################################################

#eclass

###############################################################################

class cluster_null_score_gentle:

  name   = 'ClusterNullScoreGentle';
  scores = None;
  perm   = {};
  nperm  = {};
  factor = None;

  #############################################################################

  def __init__(self, **kwargs):
    self.scores = kwargs.pop('scores');
    self.nperm  = kwargs.pop('nperm')  if 'nperm'  in kwargs else 20000;
    self.factor = kwargs.pop('factor') if 'factor' in kwargs else 1;
    perm        = {};
  #edef

  #############################################################################

  def permute(self, n, nue):
    k = (n, nue);

    if not(k in self.perm):
      S = [ 0 ] * self.nperm;
      for i in xrange(self.nperm):
        s = random.sample(self.scores, n);
        S[i] = 2*sum(s[0:n]) - (self.factor * nue);
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

class cluster_null_score_strict:

  name   = 'ClusterNullScoreStrict'
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
