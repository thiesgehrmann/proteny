#!/usr/bin/python


###############################################################################

bases = ['t', 'c', 'a', 'g'];
codons = [a+b+c for a in bases for b in bases for c in bases];
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
codon_table = dict(zip(codons, amino_acids));

###############################################################################

def translate(seq):
  aa = '';
  for i in xrange(0, len(seq), 3):
    cod = seq[i:i+3].lower();
    aa += codon_table[cod] if (cod in codon_table) else '*';
  #efor
  return aa;
#edef

###############################################################################

def revcomp(seq):
  dict = { 'A' : 'T', 'T' : 'A', 'a' : 't', 't' : 'a',
           'C' : 'G', 'G' : 'C', 'c' : 'g', 'g' : 'c',
           'R' : 'Y', 'Y' : 'R', 'r' : 'y', 'y' : 'r',
           'K' : 'M', 'M' : 'K', 'k' : 'm', 'm' : 'k',
           'S' : 'W', 'W' : 'S', 's' : 'w', 'w' : 's',
           'B' : 'V', 'V' : 'B', 'b' : 'v', 'v' : 'b',
           'D' : 'H', 'H' : 'D', 'd' : 'h', 'h' : 'd',
           'N' : 'N', 'X' : 'X', 'n' : 'n', 'x' : 'x' }
  
  return ''.join([ dict[b] for b in seq[::-1] ]);

#edef

###############################################################################
