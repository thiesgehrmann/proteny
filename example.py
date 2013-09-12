#!/usr/bin/python

import proteny as pp;
import smoothing;
import data;

sc = data.schco2();
ag = data.agabi();

PR = pp.proteny();

id_a = PR.add_org(*sc, isfile=False);
id_b = PR.add_org(*ag, isfile=False);

k = PR.blast(id_a=id_a, id_b=id_b);

PR.windows(k);

window = smoothing.windows(PR.blast_hits[k]);
clusts = smoothing.cluster(PR.blast_hits[k]);


