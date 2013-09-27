# Estimate a NULL distribution for significance testing.

import numpy as np;
import matplotlib.pyplot as P;
import math;

###############################################################################

def gene_distances(PR, k):
  """dists = gene_distances(PR, k)
     PR: The proteny instance
     k:  The task key

     Outputs:
       dists: A list of distances between genes
  """
  # Average distances between genes.

  dists = [];
  for i in [k['id_a'], k['id_b']]:
    E = PR.org_exons[i]
    odists = []
    R = [];
    G = E.Without(_.sequence, _.transcriptid, _.exonid, _.strand).GroupBy(_.geneid).Get(_.chrid[0], _.start.Min(), _.end.Max(), _.geneid).GroupBy(_.chrid);
    for (chrid, start, end, geneid) in zip(*G()):
      r = zip(start, end);
      r.sort(key=lambda x: x[1]);
      R += r;
      for i in xrange(len(r)-1):
        d = r[i+1][0] - r[i][1];
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

def getpercentile(dists):
  """pb_values = getpercentile(dists, p)
     dists: Output of gene_distances
     p:     Percentile

     Outputs:
       pb_values = The value at the 95th percentile of the distances between genes
  """

  pb_values = [];
  for d in dists:
    d = sorted(d);
    pb_values.append(d[int(len(d)*p-1)]);
  #efor

  return pb_values;
#edef

###############################################################################

def tc_threshold(pb_values):
  """thresh = tc_threshold(pb_values)
     pb_values: Output from getpercentiles

     Outputs:
       Sum of pb_values
  """
  return sum(pb_values);
#edef

###############################################################################

def get_tc_threshold(PR, k, p=0.95):
  return tc_threshold(getpercentile(gene_distances(PR, k), p));
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

