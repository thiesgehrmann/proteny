#!/bin/sh

###############################################################################

function usage(){

  arg0="$0";

  echo "Run the proteny software";
  echo "  Usage: $arg0 <name_org1> <genes_org1> <genome_org1> <name_org2> <genes_org2> <genome_org2> <p-value> <c-thresh> <outdir>";
  echo "";
  echo "  Arguments:";
  echo "    <name_orgX>:   The short name for organism X";
  echo "    <genes_orgX>:  The gene description file for organism X";
  echo "    <genome_orgX>: The genome fasta file for organism X";
  echo "    <p-value>:     The p-value threshold for clusters";
  echo "    <c-thresh>:    The conservation threshold for clusters";
  echo "    <outdir>:      The output directory to use";
  echo "    <ncores>:      The number of cores available to the ipcluster";

}

###############################################################################

if [ "$#" -ne 10 ]; then
  echo $@;
  usage $0;
  exit 1;
fi

###############################################################################

SCRIPTDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

name_org1="$1"
name_org2="$4"

genes_org1="$2";
genes_org2="$5";

genome_org1="$3";
genome_org2="$6";

pvalue="$7";
cthresh="$8";
outdir="$9";
ncores="${10}";

###############################################################################
echo "STARTING IPcluster";
export PYTHONPATH=$SCRIPTDIR:$PYTHONPATH;
ipcluster start -n "$ncores" &
sleep 10;

###############################################################################
echo "RUNNING PROTENY";
mkdir -p "$outdir/data"; 
time "${SCRIPTDIR}/run_proteny.py" "$name_org1" "$genes_org1" "$genome_org1" "$name_org2" "$genes_org2" "$genome_org2" "$pvalue" "$cthresh" "$outdir";

###############################################################################
echo "RUNNING CIRCOS"
cp "${SCRIPTDIR}/circos_run" "$outdir";
cd "$outdir";
./circos_run;

###############################################################################
echo "COMPLETE";
