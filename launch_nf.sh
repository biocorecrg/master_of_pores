#!/usr/bin/env bash
#SBATCH --no-requeue
#SBATCH --mem 7G
#SBATCH -p genoa64
#SBATCH --qos pipelines

set -e
set -u

_term() {
        echo "Caught SIGTERM signal!"
        kill -s SIGTERM $pid
        wait $pid
}

trap _term TERM

export NXF_JVM_ARGS="-Xms2g -Xmx6g"

"$@" & pid=$!

echo "Waiting for ${pid}"
wait $pid

exit 0
