#!/usr/bin/env python

import json
import sys
import inspect, os

import utils as utils
import gff3utils as Gutils

###############################################################################

__INSTALL_DIR__ = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

if len(sys.argv) < 3:
  print("Error, incomplete arguments to checkInput.py")
  sys.exit(1)
#fi

configFile = sys.argv[1];
action     = sys.argv[2];

errors = []
warnings = []

config = {}
if not(os.path.isfile(sys.argv[1])):
  errors.append("ConfigFile '%s'doesn't exist!"% sys.argv[1])
else:
  try:
    with open(sys.argv[1],"r") as ifd:
      config = json.load(ifd)
    #ewith
  except json.decoder.JSONDecodeError as error:
    errors.append("JSON error: {0}".format(error))
    config = {}
  except:
    errors.append("Could not read JSON file.")
    config = {}
  #etry
#fi
dconfig = json.load(open("%s/defaults.json" % __INSTALL_DIR__,"r"))

###############################################################################

if "genomes" not in config:
  errors.append("No 'genomes' entry in the JSON file. No genomes defined.")
else:
  if set(config["genomes"].keys()) != set([ "genome_1", "genome_2" ]):
    errors.append("You MUST define an entry for 'genome_1' and 'genome_2' in your JSON config file")
  else:
    usedNames = []
    for genome in config["genomes"]:
      if set(config["genomes"][genome].keys()) != set(["name", "genome", "gff"]):
        print(set(config["genomes"][genome].keys()))
        errors.append("You must define a 'name', 'genome' and 'gff' field for '%s'" % genome)
      else:
        if (config["genomes"][genome]["name"] == "") or (config["genomes"][genome]["name"] in usedNames):
          errors.append("You must define a unique, non-empty 'name' field for '%s'" % genome)
        #fi
        usedNames.append(config["genomes"][genome]["name"])

        fastaPassed = True
        gffPassed = True

        if (config["genomes"][genome]["genome"] == "") or not(os.path.isfile(config["genomes"][genome]["genome"])):
          errors.append("The given genome file '%s' for genome '%s' does not exist." % (config["genomes"][genome]["genome"], genome))
          fastaPassed = False
        #fi

        if (config["genomes"][genome]["gff"] == "") or not(os.path.isfile(config["genomes"][genome]["gff"])):
           fastaPassed = False
           errors.append("The given gff file '%s' for genome '%s' does not exist." % (config["genomes"][genome]["gff"], genome))
        #fi

        if fastaPassed and gffPassed:
          F = utils.loadFasta(config["genomes"][genome]["genome"])
          G = Gutils.readGFF3File(config["genomes"][genome]["gff"])

          if set([ s.split(' ')[0] for s in F.keys() ]) != set(list(G.seqids)):
            errors.append("The GFF file and Fasta files provided for genome '%s' do not refer to exactly the same contigs." % genome)
          #fi 
      #fi
    #efor
  #fi
#fi

###############################################################################
    
for error in errors:
  print("ERROR: %s" % error)
#efor

for warning in warnings:
  print("WARNING: %s" % warning)
#efor

if len(errors) > 0:
  sys.exit(1)
#fi

sys.exit(0)
