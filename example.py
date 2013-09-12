#!/usr/bin/python

import proteny as pp;
import data;

sc = data.schco2();
ag = data.agabi();

PR = pp.proteny();

PR.add_org(*sc, isrep=True);
PR.add_org(*ag, isrep=True);

PR.blast();


