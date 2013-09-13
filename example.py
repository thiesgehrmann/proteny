#!/usr/bin/python

import proteny as pp;
import data;

reload(pp);

sc = data.schco2();
ag = data.agabi();

reload(pp);
PR = pp.proteny();

id_a = PR.add_org(*sc, isfile=False);
id_b = PR.add_org(*ag, isfile=False);

k = PR.blast(id_a=id_a, id_b=id_b);

PR.save('testPR.pickle');

k = PR.windows(k);
PR.save('testPR.pickle');

k = PR.cluster_distance(k);
PR.save('testPR.pickle');

k = PR.cluster_linkage(k);
PR.save('testPR.pickle');

k = PR.cluster_hits(k);
PR.save('testPR.pickle');



