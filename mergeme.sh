#!/bin/bash
#Cq1B split into 8 parts
sbatch -A ibex-cs -J Cq1B -o gCq1B.out -e gCq1B.err --time=1-00:00:00 --account=c2071 --wrap="sh ./merge_all_using_GATK.sh Cq1B 8"
#Cq1A split into 7 parts
sbatch -A ibex-cs -J Cq1A -o gCq1A.out -e gCq1A.err --time=1-00:00:00 --account=c2071 --wrap="sh ./merge_all_using_GATK.sh Cq1A 7"
#Cq4A split into 7 parts
sbatch -A ibex-cs -J Cq4A -o gCq4A.out -e gCq4A.err --time=1-00:00:00 --account=c2071 --wrap="sh ./merge_all_using_GATK.sh Cq4A 7"
#Cq5B split into 9 parts
sbatch -A ibex-cs -J Cq5B -o gCq5B.out -e gCq5B.err --time=1-00:00:00 --account=c2071 --wrap="sh ./merge_all_using_GATK.sh Cq5B 9"
#Cq5A split into 8 parts
sbatch -A ibex-cs -J Cq5A -o gCq5A.out -e gCq5A.err --time=1-00:00:00 --account=c2071 --wrap="sh ./merge_all_using_GATK.sh Cq5A 8"
#Cq6A split into 8 parts
sbatch -A ibex-cs -J Cq6A -o gCq6A.out -e gCq6A.err --time=1-00:00:00 --account=c2071 --wrap="sh ./merge_all_using_GATK.sh Cq6A 8"
#Cq7A split into 7 parts
sbatch -A ibex-cs -J Cq7A -o gCq7A.out -e gCq7A.err --time=1-00:00:00 --account=c2071 --wrap="sh ./merge_all_using_GATK.sh Cq7A 7"
#Cq2A split into 6 parts
sbatch -A ibex-cs -J Cq2A -o gCq2A.out -e gCq2A.err --time=1-00:00:00 --account=c2071 --wrap="sh ./merge_all_using_GATK.sh Cq2A 6"
#Cq2B split into 9 parts
sbatch -A ibex-cs -J Cq2B -o gCq2B.out -e gCq2B.err --time=1-00:00:00 --account=c2071 --wrap="sh ./merge_all_using_GATK.sh Cq2B 9"
#Cq7B split into 9 parts
sbatch -A ibex-cs -J Cq7B -o gCq7B.out -e gCq7B.err --time=1-00:00:00 --account=c2071 --wrap="sh ./merge_all_using_GATK.sh Cq7B 9"
#Cq4B split into 8 parts
sbatch -A ibex-cs -J Cq4B -o gCq4B.out -e gCq4B.err --time=1-00:00:00 --account=c2071 --wrap="sh ./merge_all_using_GATK.sh Cq4B 8"
#Cq3A split into 7 parts
sbatch -A ibex-cs -J Cq3A -o gCq3A.out -e gCq3A.err --time=1-00:00:00 --account=c2071 --wrap="sh ./merge_all_using_GATK.sh Cq3A 7"
#Cq9A split into 7 parts
sbatch -A ibex-cs -J Cq9A -o gCq9A.out -e gCq9A.err --time=1-00:00:00 --account=c2071 --wrap="sh ./merge_all_using_GATK.sh Cq9A 7"
#Cq3B split into 9 parts
sbatch -A ibex-cs -J Cq3B -o gCq3B.out -e gCq3B.err --time=1-00:00:00 --account=c2071 --wrap="sh ./merge_all_using_GATK.sh Cq3B 9"
#Cq9B split into 8 parts
sbatch -A ibex-cs -J Cq9B -o gCq9B.out -e gCq9B.err --time=1-00:00:00 --account=c2071 --wrap="sh ./merge_all_using_GATK.sh Cq9B 8"
#Cq8B split into 9 parts
sbatch -A ibex-cs -J Cq8B -o gCq8B.out -e gCq8B.err --time=1-00:00:00 --account=c2071 --wrap="sh ./merge_all_using_GATK.sh Cq8B 9"
#Cq6B split into 10 parts
sbatch -A ibex-cs -J Cq6B -o gCq6B.out -e gCq6B.err --time=1-00:00:00 --account=c2071 --wrap="sh ./merge_all_using_GATK.sh Cq6B 10"
#Cq8A split into 7 parts
sbatch -A ibex-cs -J Cq8A -o gCq8A.out -e gCq8A.err --time=1-00:00:00 --account=c2071 --wrap="sh ./merge_all_using_GATK.sh Cq8A 7"
#CqUO split into 15 parts
sbatch -A ibex-cs -J CqUO -o gCqUO.out -e gCqUO.err --time=1-00:00:00 --account=c2071 --wrap="sh ./merge_all_using_GATK.sh CqUO 15"
#MK159176 split into 1 part
sbatch -A ibex-cs -J MK159176 -o gMK159176.out -e gMK159176.err --time=1-00:00:00 --account=c2071 --wrap="sh ./merge_all_using_GATK.sh MK159176 1"
#MK182703 split into 1 part
sbatch -A ibex-cs -J MK182703 -o gMK182703.out -e gMK182703.err --time=1-00:00:00 --account=c2071 --wrap="sh ./merge_all_using_GATK.sh MK182703 1"
