import sys;
sys.path.append('utils');
sys.path.append('visualization');

import cluster_null as null;
import circos_data as cdata;
import circos_region as cr;
import circos_chr    as cc;
import proteny as ps;
import data;

###############################################################################

savename = 'PR_paths.proteny';
circdir  = 'circos';

###############################################################################

pleos = data.pleos2();
schco = data.schco2();
agabi = data.agabi2();


a = PR.add_org(*schco, isfile=False);
b = PR.add_org(*pleos, isfile=False);
c = PR.add_org(*agabi, isfile=False);

k_ab = PR.analyze(id_a=a, id_b=b);
k_bc = PR.analyze(id_a=b, id_b=c);
k_ac = PR.analyze(id_a=a, id_b=c);

PR.save(savename);


