#!/bin/bash


export SCRIPT=/ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/downstream_analysis_generic_v2.0_ALL_allsites.sh
export REF=/ibex/scratch/projects/c2071/1000quinoa/reference/quinoa_pb_chicago-2-final_PBJELLY_pilon_pseudo_RENAMED_CqUO_Cp_Mt.fasta
export INPUT=/ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/VCF_MERGE
export OUTPUT=/ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/OUTPUT
$SCRIPT $REF $INPUT $OUTPUT

