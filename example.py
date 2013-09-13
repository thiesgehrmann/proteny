#!/usr/bin/python

import proteny as pp;
import data;

savename = 'PR.proteny';

sc = data.schco2();
ag = data.agabi();

PR = pp.proteny();

id_a = PR.add_org(*sc, isfile=False);
id_b = PR.add_org(*ag, isfile=False);

k = PR.blast(id_a=id_a, id_b=id_b);

PR.save(savename);

k = PR.windows(k);
PR.save(savename);

k = PR.hit_distance(k);
PR.save(savename);

k = PR.cluster_linkage(k);
PR.save(savename);

k = PR.cluster_hits(k);
PR.save(savename);



