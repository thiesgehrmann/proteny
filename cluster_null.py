from scipy import stats;
import random;
import bisect;
import numpy as np;
import math;
import copy;
import time;
import util;
import types;
import cPickle as cpickle;

from IPython.parallel import Client

###############################################################################

#class cluster_null_clt:
#
#  name   = 'ClusterNullCLT';
#  H_mu   = None;
#  H_s2   = None;
#  U1_mu  = None;
#  U1_s2  = None;
#  U2_mu  = None;
#  U2_s2  = None;
#
#  #############################################################################
#
#  def __init__(self, **kwargs):
#    scores = kwargs.pop('scores');
#    hits   = kwargs.pop('hits');
#    hitA   = kwargs.pop('hitA');
#    hitB   = kwargs.pop('hitB');    
#
#    self.H_mu  = self.estmu(scores);
#    self.H_s2  = self.ests2(scores);
#
#    L          = [ hits[hitA[k]][12] for k in hitA.keys() ];
#    self.U1_mu = self.estmu(L);
#    self.U1_s2 = self.ests2(L);
#
#    L          = [ hits[hitB[k]][12] for k in hitB.keys() ];
#    self.U2_mu = self.estmu(L);
#    self.U2_s2 = self.ests2(L);
#  #edef
#
#  #############################################################################
#
#  def estmu(self, L):
#    mu = np.mean(L);
#    return mu;
#  #edef
#
#  #############################################################################
#
#  def ests2(self, L):
#    s2 = np.var(L);
#    return s2;
#  #edef
#
#  #############################################################################
#
#  def pvalue(self, **kwargs):
#    # The score is a sum of CLT estimated normal distributions.
#
#    score = kwargs.pop('score');
#    n     = kwargs.pop('n');
#    nue_a = kwargs.pop('nue_a');
#    nue_b = kwargs.pop('nue_b');
#
#    mu = (self.H_mu * n) - (self.U1_mu * nue_a) - (self.U2_mu * nue_b);
#    s2 = (self.H_s2 * n) + (self.U1_s2 * nue_a) + (self.U2_s2 * nue_b);
#
#    z = ( score - mu ) / s2;
#    p = 1 - stats.norm.cdf(z);
#
#    return p;
#  #edef
#
#  #############################################################################
#
##eclass

###############################################################################

#class cluster_null_score_gentle:
#
#  name   = 'ClusterNullScoreGentle';
#  scores = None;
#  perm   = {};
#  nperm  = {};
#
#  #############################################################################
#
#  def __init__(self, **kwargs):
#    self.scores = kwargs.pop('scores');
#    self.nperm  = kwargs.pop('nperm')  if 'nperm'  in kwargs else 20000;
#    perm        = {};
#  #edef
#
#  #############################################################################
#
#  def permute(self, n, nue):
#    k = (n, nue);
#
#    if not(k in self.perm):
#      S = [ 0 ] * self.nperm;
#      for i in xrange(self.nperm):
#        s = random.sample(self.scores, n);
#        S[i] = 2*sum(s[0:n]) - (nue);
#      #efor
#      self.perm[k] = sorted(S);
#    #fi
#
#    return self.perm[k];
#  #edef
#
#  #############################################################################
#
#  def pvalue(self, **kwargs):
#    score = kwargs.pop('score');
#    n     = kwargs.pop('n');
#    nue_a = kwargs.pop('nue_a');
#    nue_b = kwargs.pop('nue_b');
#
#    nue = nue_a + nue_b;
#
#    if nue >= n:
#      return 1.0;
#    #fi
#
#    p = self.permute(n, nue);
#
#    i = bisect.bisect_left(p, score);
#
#    return 1.0 - float(i) / float(self.nperm);
#  #edef
#
#  #############################################################################
#
##eclass
#
################################################################################
#
#class cluster_null_score_strict:
#
#  name   = 'ClusterNullScoreStrict'
#  scores = None;
#  hits   = None;
#  hitA   = None;
#  hitB   = None;
#  perm   = {};
#  nperm  = {};
#
#  #############################################################################
#
#  def __init__(self, **kwargs):
#    self.scores = kwargs.pop('scores');
#    self.nperm  = kwargs.pop('nperm') if 'nperm' in kwargs else 10000000;
#    self.hits   = kwargs.pop('hits');
#    self.hitA   = kwargs.pop('hitA');
#    self.hitB   = kwargs.pop('hitB');
#    perm        = {};
#  #edef
#
#  #############################################################################
#
#  def permute(self, n, nue_a, nue_b):
#    k = (n, nue_a, nue_b);
#    
#    if not(k in self.perm):
#      S = [ 0 ] * self.nperm;
#      for i in xrange(self.nperm):
#        s = random.sample(self.hits, n + nue_a + nue_b);
#        match_s = 2 * sum([x[12] for x in s[0:n]]);
#        ue_a_s  = sum([self.hits[self.hitA[(x[2], x[3])]][12] for x in s[n:n+nue_a]])
#        ue_b_s  = sum([self.hits[self.hitB[(x[5], x[6])]][12] for x in s[n+nue_a:]])
#        S[i] = match_s - ue_a_s - ue_b_s;
#      #efor
#      self.perm[k] = sorted(S);
#    #fi
#    
#    return self.perm[k];
#  #edef
#
#  #############################################################################
#
#  def pvalue(self, **kwargs):
#    score = kwargs.pop('score');
#    n     = kwargs.pop('n');
#    nue_a = kwargs.pop('nue_a');
#    nue_b = kwargs.pop('nue_b');
#
#    if nue_a + nue_b >= n:
#      return 1.0;
#    #fi
#    
#    p = self.permute(n, nue_a, nue_b);
#    
#    i = bisect.bisect_left(p, score) + 1;
#    
#    return 1.0 - float(i) / float(self.nperm + 1);
#  #edef
#
#  #############################################################################
#
##eclass

###############################################################################


class cluster_null_score_strict_smart:

  # With idea from:
  # Knijnenburg, T. a, Wessels, L. F. a, Reinders, M. J. T., & Shmulevich, I. (2009). Fewer permutations, more accurate P-values. Bioinformatics (Oxford, England)

  name   = 'ClusterNullScoreStrictSmart'
  scores = None;
  hits   = None;
  hitA   = None;
  hitB   = None;
  perm   = {};
  sperm  = None; # Steps of permutations necessary
  mperm  = None;
  perm_func = None;
  parc      = None;
  nproc     = 8;

  M           = None;
  Mex         = 10;
  psuedocount = 1;

  scores = None;
  hitAs  = None;
  hitBs  = None;

  clt_H_mu   = None;
  clt_H_s2   = None;
  clt_U1_mu  = None;
  clt_U1_s2  = None;
  clt_U2_mu  = None;
  clt_U2_s2  = None;

  #############################################################################

  def __init__(self, **kwargs):

    if 'storage' in kwargs and util.fex(kwargs['storage']):
      self.init_from_file(**kwargs);
    else:
      self.init_now(**kwargs);
    #fi
  #edef

  #############################################################################

  def init_from_file(self, **kwargs):
    D = cpickle.load(open(kwargs['storage'], 'r'));
    self.perm = D;
    self.init_now(**kwargs);
  #edef

  #############################################################################
    
  def store_perms(self):
    f = open(self.storage, 'w');
    cpickle.dump(self.perm, f);
    f.close();
  #edef

  #############################################################################

  def init_now(self, **kwargs):

    self.storage = kwargs.pop('storage') if 'storage' in kwargs else 'cluster_null_score_strict_smart_storage.dat';

    self.sperm  = kwargs.pop('sperm') if 'sperm' in kwargs else 1000;
    self.M      = self.sperm;
    self.alpha  = kwargs.pop('alpha');
    self.tests  = kwargs.pop('tests');
    self.mperm  = math.ceil(1.0/(self.alpha / self.tests));

    hits   = kwargs.pop('hits');
    hitA   = kwargs.pop('hitA');
    hitB   = kwargs.pop('hitB');

    self.nproc     = kwargs.pop('nproc') if 'nproc' in kwargs else 8;
    self.parthresh = kwargs.pop('parthresh') if 'parthresh' in kwargs else 25000;

    self.psuedocount = 0 if ('psuedocount' in kwargs and kwargs['psuedocount'] == False) else 1;

    self.scores = np.array([ h[12] for h in hits ], dtype=float);
    self.hitAs  = np.array([hits[v][12] for v in hitA.values()], dtype=float);
    self.hitBs  = np.array([hits[v][12] for v in hitB.values()], dtype=float);

    self.clt_H_mu   = np.mean(self.scores);
    self.clt_H_s2   = np.var(self.scores);
    self.clt_U1_mu  = np.mean(self.hitAs);
    self.clt_U1_s2  = np.var(self.hitAs);
    self.clt_U2_mu  = np.mean(self.hitBs);
    self.clt_U2_s2  = np.var(self.hitBs);

    p = util.run_cmd("ipcluster start -n %d" % (self.nproc), bg=True)

    i = 0;
    self.perm_func = self.permute;
    for i in xrange(120):
      try:
        rc             = Client();
        self.parc      = rc.load_balanced_view();
        self.perm_func = self.permute_par;
        print 'Parallel option';
        break;
      except IOError:
        time.sleep(1)
        pass;
      #etry
    #ewhile

  #edef

  #############################################################################

  def update_mperm(self, tests):
    self.tests = tests
    self.mperm = math.ceil(1.0/(self.alpha / self.tests));
  #edef

  #############################################################################

  def permute(self, n, nue_a, nue_b, more=False):
    k = (n, nue_a, nue_b);

    if not(k in self.perm):
      self.perm[k] = (0, []);
      more = True;
    else:
      print "HIT!";
    #fi

    perms_done, topM = self.perm[k];

    if more:
      S = [ 0 ] * self.sperm;
      indx = np.arange(len(self.hits), dtype=int);
      for i in xrange(self.sperm):
        s_indx = np.array(random.sample(indx, (n + nue_a + nue_b)), dtype=int);
                                                  
        match_s = 2 * self.scores[s_indx[0:n]].sum();
        ue_a_s  = self.hitAs[s_indx[n:nue_a]].sum();
        ue_b_s  = self.hitBs[s_indx[n+nue_a:]].sum();

        S[i] = match_s - ue_a_s - ue_b_s;
      #efor
      total_perms_done = perms_done + self.sperm;
      new_topM         = sorted(S + topM)[-self.M:];
      self.perm[k] = (total_perms_done, new_topM);
    #fi
    
    return self.perm[k];
  #edef

  #############################################################################

  def permute_par(self, n, nue_a, nue_b, more=False):
    k = (n, nue_a, nue_b);

    if not(k in self.perm):
      self.perm[k] = (0, []);
      more = True;
    else:
      print "HIT";
    #fi
    perms_done, topM = self.perm[k];

    nperm = max(perms_done, self.sperm);

    if (2*nperm > self.mperm) and (self.mperm > self.sperm):
      nperm = min(nperm, abs(self.mperm - nperm + 1));
    #fi

    J = [None] * self.nproc;
    P = [];
    if more:
      if nperm > self.parthresh:
        print "\nPerforming %d permutations" % nperm;
        P = permute_par_outer(self.parc, n, nue_a, nue_b, nperm, self.nproc, self.scores, self.hitAs, self.hitBs, self.M);
      else:
        P = permute_par_core(n, nue_a, nue_b, int(nperm), self.scores, self.hitAs, self.hitBs, self.M);
      #fi
      total_perms_done = perms_done + nperm;
      sort_perms       = sorted(P + topM);
      new_topM         = sort_perms[-self.M:]
      self.perm[k]     = (total_perms_done, new_topM);
    #fi
    return self.perm[k];
  #edef

  #############################################################################

  def pvalue(self, **kwargs):
    score = kwargs.pop('score');

    n     = kwargs.pop('n');
    nue_a = kwargs.pop('nue_a');
    nue_b = kwargs.pop('nue_b');

    if n > 10 and nue_a > 10 and nue_b > 10:
      return self.pvalue_clt(score, n, nue_a, nue_b);
    else:
      return self.pvalue_perm(score, n, nue_a, nue_b);
    #fi
  #edef

  #############################################################################

  def pvalue_perm(self, score, n, nue_a, nue_b):

    perms_done, topM = self.perm_func(n, nue_a, nue_b);

    while True:
      i = bisect.bisect_left(topM, score);
        # We have more than Mex exceedences, we don't need a psuedocount, and we can stop permuting
      if i < (self.M - self.Mex):
        pval = 1.0 - (float(perms_done - self.M + i) / float(perms_done))
        break;
        # If, however, we don't have at least Mex exceedences, we should only stop permuting if we haven't done enough permutations...
      elif perms_done >= self.mperm:
        #print "n: %d, nue_a: %d, nue_b: %d" % (n, nue_a, nue_b);
        #print "Perms Done: %d\n" % perms_done
        pval = 1.0 - (float(perms_done - self.M + i) / float(perms_done + self.psuedocount));
        break;
      else:
        perms_done, topM = self.perm_func(n, nue_a, nue_b, more=True);
      #fi
    #ewhile

    #print topM, score, i, (perms_done - self.M + i), pval;

    return (stats.norm.ppf(1.0 - pval), pval);
  #edef

  #############################################################################

  def pvalue_clt(self, score, n, nue_a, nue_b):

    mu = 2*(self.clt_H_mu * n) - (self.clt_U1_mu * nue_a) - (self.clt_U2_mu * nue_b);
    s2 = (4 * self.clt_H_s2 * n) + (self.clt_U1_s2 * nue_a) + (self.clt_U2_s2 * nue_b);

    z = ( score - mu ) / math.sqrt(s2);
    p = 1 - stats.norm.cdf(z);

    #print ((mu, s2, z), n, nue_a, nue_b, p);
    #print "MU:" + str((2*(self.clt_H_mu * n), (self.clt_U1_mu * nue_a), (self.clt_U2_mu * nue_b)));
    #print "VAR:" + str(((4 * self.clt_H_s2 * n), (self.clt_U1_s2 * nue_a), (self.clt_U2_s2 * nue_b)))

    return (z, p);
  #edef

  #############################################################################

#eclass

###############################################################################

##class ClusterNullsimple:
#
#  name   = 'ClusterNullSimple'
#  scores = None;
#  hits   = None;
#  hitA   = None;
#  hitB   = None;
#  perm   = {};
#  sperm  = None; # Steps of permutations necessary
#  mperm  = None;
#  perm_func = None;
#  parc      = None;
#  nproc     = 8;
#
#  M           = 10;
#  psuedocount = 1;
#
#  scores = None;
#  hitAs  = None;
#  hitBs  = None;
#
#  clt_H_mu   = None;
#  clt_H_s2   = None;
#  clt_U1_mu  = None;
#  clt_U1_s2  = None;
#  clt_U2_mu  = None;
#  clt_U2_s2  = None;
#
#  #############################################################################
#
#  def __init__(self, **kwargs):
#    
#    if 'storage' in kwargs and util.fex(kwargs['storage']):
#      self.init_from_file(**kwargs);
#    else:
#      self.init_now(**kwargs);
#    #fi
#  #edef
#
#  #############################################################################
#
#  def init_from_file(self, **kwargs):
#    D = pickle.load(open(kwargs['storage'], 'r'));
#    self.perm = D;
#    self.init_now(**kwargs);
#  #edef
#
#  #############################################################################
#
#  def store_perms(self):
#    f = open(self.storage, 'w');
#    pickle.dump(self.perm, f);
#    f.close();
#  #edef
#
#  #############################################################################
#
#  def init_now(self, **kwargs):
#    
#    self.storage = kwargs.pop('storage') if 'storage' in kwargs else 'cluster_null_simple_storage.dat';
#
#    self.scores = kwargs.pop('scores');
#    self.hits   = kwargs.pop('hits');
#    self.hitA   = kwargs.pop('hitA');
#    self.hitB   = kwargs.pop('hitB');
#
#    self.nperm       = kwargs.pop('nperm') if 'nperm' in kwargs else 10000;
#    self.nproc       = kwargs.pop('nproc') if 'nproc' in kwargs else 8;    
#    self.psuedocount = 0 if ('psuedocount' in kwargs and kwargs['psuedocount'] == False) else 1;
#    
#    self.scores = np.array([ h[12] for h in self.hits ], dtype=float);
#    self.hitAs  = np.array([ self.hits[self.hitA[(h[2], h[3])]][12] for h in self.hits], dtype=float);
#    self.hitBs  = np.array([ self.hits[self.hitB[(h[5], h[6])]][12] for h in self.hits], dtype=float);
#   
#    vals_hitAs  = [hits[v][12] for v in hitA.values()];
#    vals_hitBs  = [hits[v][12] for v in hitB.values()];
# 
#    self.clt_H_mu   = np.mean(self.scores);
#    self.clt_H_s2   = np.var(self.scores);
#    self.clt_U1_mu  = np.mean(vals_hitAs);
#    self.clt_U1_s2  = np.var(vals_hitAs);
#    self.clt_U2_mu  = np.mean(vals_hitBs);
#    self.clt_U2_s2  = np.var(vals_hitBs);
#    
#    p = util.run_cmd("ipcluster start -n %d" % (self.nproc), bg=True)
#    
#    i = 0;
#    for i in xrange(120):
#      try:
#        rc             = Client();
#        self.parc      = rc.load_balanced_view();
#        print 'Parallel option';
#        break;
#      except IOError:
#        time.sleep(1)
#        pass;
#      #etry
#    #ewhile
#
#  #edef
#
#  #############################################################################
#
#  def pvalue_perm(self, score, n, nue_a, nue_b):
#
#    k = (n, nue_a, nue_b);
#
#    if k in self.perm:
#      P = self.perm[k];
#    else:
#      print "Running par core now, for %d permutations" % self.nperm;
#      print "n: %d, nue_a: %d, nue_b: %d" % (n, nue_a, nue_b);
#      P = permute_par_core(n, nue_a, nue_b, self.nperm, self.scores, self.hitAs, self.hitBs)
#      self.perm[k] = P;
#    #fi
#
#    pvalue = 1.0 - (float(bisect.bisect_left(P, score)) / float(len(P)+self.psuedocount));
#
#    return pvalue;
#
#  #edef
#
#  #############################################################################
#
#  def pvalue_clt(self, score, n, nue_a, nue_b):
#    mu = (self.clt_H_mu * n) - (self.clt_U1_mu * nue_a) - (self.clt_U2_mu * nue_b);
#    s2 = (self.clt_H_s2 * n) + (self.clt_U1_s2 * nue_a) + (self.clt_U2_s2 * nue_b);
#    
#    z = ( score - mu ) / math.sqrt(s2);
#    p = 1 - stats.norm.cdf(z);
#    print ((mu, s2, z), n, nue_a, nue_b, p);
#
#    return p;
#  #edef
#
##eclass

###############################################################################

def permute_par_outer(parc, n, nue_a, nue_b, nperm, nproc, scores, hitAs, hitBs, topM):

  J = [None]*nproc;
  P = [];

  for i in xrange(nproc):
    J[i] = parc.apply_async(permute_par_core, n, nue_a, nue_b, int(math.ceil(nperm/nproc)), scores, hitAs, hitBs, topM)
  #efor
  for j in J:
    P = P + j.get();
  #efor

  print 2*n*np.mean(scores), nue_a * np.mean(hitAs), nue_b * np.mean(hitBs);
  

  return sorted(P)[-topM:];
#edef

###############################################################################

def permute_par_core(n, nue_a, nue_b, nperm, s1, s2, s3, topM):
  import numpy as np;
  import random;

  choice = np.random.choice;

  S = [ 0 ] * nperm;

  s1_indx = np.arange(len(s1), dtype=int);
  s2_indx = np.arange(len(s2), dtype=int);
  s3_indx = np.arange(len(s3), dtype=int);

  nperm_left = nperm;
  total_vals = n+nue_a+nue_b;

  print "--PERMUTING"

  j = 0;
  while nperm_left > 0:
 
    nperms_this_time = min(1000, nperm_left);
    r_s1_indx = choice(s1_indx, nperms_this_time * n);
    r_s2_indx = choice(s2_indx, nperms_this_time * nue_a);
    r_s3_indx = choice(s3_indx, nperms_this_time * nue_b);
    for i in xrange(nperms_this_time):

      match_s = s1[r_s1_indx[0:n]].sum();
      ue_a_s  = s2[r_s2_indx[0:nue_a]].sum();
      ue_b_s  = s3[r_s3_indx[0:nue_b]].sum();

      r_s1_indx = r_s1_indx[n:];
      r_s2_indx = r_s2_indx[nue_a:];
      r_s3_indx = r_s3_indx[nue_b:];

      S[j+i] = 2 * match_s - ue_a_s - ue_b_s;
    #efor
    j = j + nperms_this_time;
    nperm_left = nperm_left - nperms_this_time;
  #efor

  print "--RETURNING"

  return sorted(S)[-topM:];
#edef

###############################################################################

