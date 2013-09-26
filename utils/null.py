# Estimate a NULL distribution for significance testing.

import numpy as np;
import matplotlib.pyplot as P;
from scipy.stats import poisson;
import math;

###############################################################################

def gene_distances(PR):
  # Average distances between genes.

  dists = [];
  for (i, E) in enumerate(PR.org_exons):
    odists = []
    R = [];
    G = E.Without(_.sequence, _.transcriptid, _.exonid, _.strand).GroupBy(_.geneid).Get(_.chrid[0], _.start.Min(), _.end.Max(), _.geneid).GroupBy(_.chrid);
    for (chrid, start, end, geneid) in zip(*G()):
      r = zip(start, end);
      r.sort(key=lambda x: x[1]);
      R += r;
      for i in xrange(len(r)-1):
        d = r[i+1][0] - r[i][1];
        if d > 20000:
          print r[i+1], r[i];
        #fi
        if d > 0:
          odists.append(r[i+1][0] - r[i][1]);
        #fi
      #efor
    #efor
    dists.append(odists);
  #efor

  return dists;
#edef

###############################################################################

def getpercentile(dists, p=0.95):
  pb_values = [];
  for d in dists:
    d = sorted(d);
    pb_values.append(d[int(len(d)*p-1)]);
  #efor

  return pb_values;
#edef

###############################################################################

def tc_threshold(pb_values):
  return sum(pb_values);
#edef

###############################################################################

def visualize_dists(dists, cumul=False):
  #mus  = [ np.mean(d) for d in dists ];
  #stds = [ np.std(d) for d in dists ];

  if len(dists) < 4:
    x = len(dists);
    y = 1;
  else:
    x = math.ceil(sqrt(len(dists)));
    y = x;
  #fi

  for i in xrange(len(dists)):
    P.subplot(x,y, i+1);
    if cumul == True:
      n, bins, patches = P.hist(dists[i], bins=500, cumulative=True, histtype='stepfilled');
    else:
      n, bins, patches = P.hist(dists[i], bins=500, histtype='stepfilled');
    #fi
    P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
  #efor
  P.show();
#edef

###############################################################################

