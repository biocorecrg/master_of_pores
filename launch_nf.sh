#!/bin/bash
#SBATCH --no-requeue
#SBATCH --mem 8192M
#SBATCH -p genoa64
#SBATCH --qos='pipelines'
set -x
pid=""
	kill_func() {
	echo TRAP; 
	kill $pid ; 
	wait $pid
}
trap kill_func INT
trap kill_func EXIT

"$@" & pid=$! ; echo "waiting for ${pid}" ; wait $pid 
