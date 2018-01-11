#!/bin/bash
#SBATCH -J "anneal"        	 	#  Give the job a name
#SBATCH -n 1				# use 1 node
#SBATCH -d 32				# use 32 cpus
#SBATCH -t 99:00:00             	# time needed 
#SBATCH --mail-type=NONE                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sebastian.menze@uib.no   # email to user
#SBATCH --output=parallel_annealing.out     #Standard output and error log

cd /work/users/sme048/weddell_sea_scenarios
module load python
aprun -n 1 -d 32 python paraanneal.py