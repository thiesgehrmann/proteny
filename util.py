from ibidas import *;

import numpy as np;
import bisect;
import shlex, subprocess, os;
import os.path;

###############################################################################

def chr_pair_group(hits):
  """HC = chr_pair_group(H)
     H: Output list from sort_org

     Group hits by chromosomes they exist on
     NOTE: DOES NOT PRESERVE INDICES OF O LISTS!
     Outputs:
       HC: A dictionary of hits per pair of chromosomes
  """

  H_chrs = {};

  for (i, hit) in enumerate(hits):
    a_chr1 = hit[1];
    a_chr2 = hit[4];
    k = (a_chr1, a_chr2);
    if k not in H_chrs:
      H_chrs[k] = [];
    #fi
    H_chrs[k].append(i);
  #efor

  return H_chrs;

#edef

###############################################################################

def prep_exon_list(E):
  L = zip(*E.Get(_.chrid, _.start, _.transcriptid, _.exonid)());

  chrs   = {};
  for (chrid, start, transcriptid, exonid) in L:
    if chrid not in chrs:
      chrs[chrid]   = [];
    #fi
    chrs[chrid].append((start,transcriptid,exonid));
  #efor

  for chrid in chrs.keys():
    chrs[chrid].sort(key=lambda x: x[0]);
  #efor

  return chrs;
#edef

###############################################################################

def exons_in_reg(chrs, chrid, start, end):
  E = chrs[chrid];
  L = [ x[0] for x in E ];
  i = bisect.bisect_right(L, start);
  j = bisect.bisect_left(L, end);

  return [ (k[1], k[2]) for k in E[i-1:j] ];

#edef

###############################################################################

def reindex_blast(br):
  """ B + reindex_blast(br)
      br: A rep of BLAST results.
          Must contain at least:
          _.a_chrid, _.a_geneid, _.a_exonid
          _.b_chrid, _.b_geneid, _.b_exonid
          _.a_start, _.a_end
          _.b_start, _.b_end
          _.coverage, _.pident, _.evalue, _.bitscore

      Output:
        B: A Rep of blast results with an index
  """

  nbr = br.Shape()();
  ind = [i for i in xrange(nbr)];

  F = br.Get(_.a_chrid, _.a_geneid, _.a_exonid, \
             _.b_chrid, _.b_geneid, _.b_exonid, \
             _.a_start,                         \
             _.a_end,                           \
             _.b_start,                         \
             _.b_end,                           \
             _.coverage,                        \
             _.score,                           \
             _.pident,                          \
             _.evalue,                          \
             _.bitscore);
  F = Rep(zip(ind, *F()))
  F = F / ('i', 'a_chrid', 'a_geneid', 'a_exonid', 'b_chrid', 'b_geneid', 'b_exonid', 'a_start', 'a_end', 'b_start', 'b_end', 'coverage', 'score', 'pident', 'evalue', 'bitscore');

  return F.Copy();
#edef

###############################################################################

def overlap(r1, r2):
  """
    d = overlap(r1, r2)

      Returns the overlap d of two regions specified.
      If the regions do not overlap, then the negative distance between them is returned.

      r1 = (start1, end1)
      r2 = (start2, end2);

    s1      e1
    |       |
    ####1####
           ||  <---------  (overlap)
           ####2####
           |       |
           s2      e2
  """

  if r1[0] < r2[0]:
    e1 = r1[1];
    s2 = r2[0];
  elif r1[0] == r2[0]:
    e1 = max(r1[1], r2[1]);
    s2 = r1[0];
  else:
    e1 = r2[1];
    s2 = r1[0];
  #fi

  ov  = e1 - s2;
  mov = min(r1[1] - r1[0], r2[1] - r2[0]);

  return min(mov, ov);
#edef

###############################################################################

def shuffle_sub(L, s=0, e=-1):
  NL = [];
  n  = len(L);

  index = [ i for i in xrange(n) ];
  shfld = [ i for i in xrange(n) ];
  np.random.shuffle(shfld);

  for (i, j) in zip(index, shfld):
    k = L[i][0:s] + L[j][s:e] + L[i][e:];
    NL.append(k);
  #efor

  return NL;
#edef

###############################################################################

def run_par_cmds(cmd_list, max_threads=12, stdin=None, stdout=None, stderr=None):

  p = [];
  i = 0;
  retval = 0;
  cmds = len(cmd_list);

  while True:
    while len(p) < max_threads and i < cmds:
      print "RUNNING: %s" % cmd_list[i]; sys.stdout.flush();
      p.append( (run_cmd(cmd_list[i], bg=True, stdin=stdin, stdout=stdout, stderr=stderr),i) );
      i = i + 1;
    #ewhile

    time.sleep(0.5);

    running = [ (j, k) for (j,k) in p if j.poll() == None ];
    completed = [ (j, k) for (j,k) in p if j.poll() != None ];

    for (j,k) in completed:
      if j.returncode != 0:
        retval = retval + j.returncode;
        print "ERROR: Failed in cmd: %s" % cmd_list[k]; sys.stdout.flush();
      else:
        print "COMPLETED: cmd : %s" % cmd_list[k]; sys.stdout.flush();
      #fi
    #efor
    p = running;
    if len(p) == 0:
      break;
    #fi
  #ewhile

  return retval;
#edef


###############################################################################

def run_seq_cmds(cmd_list, stdin=None, stdout=None, stderr=None):

  for cmd in [ x for x in cmd_list if x ]:
    retval = run_cmd(cmd, stdin=stdin, stdout=stdout, stderr=stderr);
    if retval != 0:
      print "ERROR: Failed on cmd: %s" % cmd;
      return retval;
    #fi
    print "COMPLETED: cmd : %s" % cmd;
    sys.stdout.flush();
  #efor

  return 0;
#edef

###############################################################################

def run_cmd(cmd, bg=False, stdin=None, stdout=None, stderr=None):
  p = subprocess.Popen(shlex.split(cmd), stdin=stdin, stdout=stdout, stderr=stderr);
  if bg:
    return p;
  else:
    (pid, r) = os.waitpid(p.pid, 0);
    return r;
  #fi
#edef

###############################################################################

def fex(fname):
  return os.path.isfile(fname);
#edef

###############################################################################

