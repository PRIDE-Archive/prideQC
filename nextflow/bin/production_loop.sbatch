#!/usr/bin/env bash
#SBATCH --job-name=prideqc-loop
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=7-00:00:00
# Codon: run production_loop.py as ONE Slurm job; the Nextflow controller of each batch
# runs on this compute node (verified 2026-11-07: Slurm submit from a compute node works).
# Usage (all loop arguments are passed through; PRIDEQC_SOURCE is the pinned source dir):
#   sbatch --export=ALL,PRIDEQC_SOURCE=<src> --output=<loopdir>/slurm-%j.log \
#     <src>/nextflow/bin/production_loop.sbatch --root ... --loop-id ... --source <src> --batches N ...
# Stop after the current batch: touch <root>/work/loops/<loop-id>/STOP
# No set -e / pipefail. The exit code is the loop result (0 DONE, 3 STOP, 2 BLOCKED).
JDK=${PRIDEQC_JDK:-/hps/software/users/juan/pride/java/corretto-17.0.6}
export JAVA_HOME=$JDK JAVA_CMD=$JDK/bin/java NXF_JAVA_HOME=$JDK NXF_VER=${NXF_VER:-26.10.2}
export PATH="$JDK/bin:${PRIDEQC_NXF_BIN_DIR:-$HOME/.local/bin}:/ebi/slurm/codon/scripts/user:/ebi/slurm/codon/bin:/usr/bin:/bin:/usr/sbin:/sbin"
umask 027
echo "LOOP_HOST=$(hostname -f) JOB=${SLURM_JOB_ID:-none} DATE=$(date -u +%FT%TZ)"
[ -f "${PRIDEQC_SOURCE:-}/nextflow/bin/production_loop.py" ] || { echo "LOOP=BLOCKED PRIDEQC_SOURCE unset or wrong"; exit 2; }
python3 "$PRIDEQC_SOURCE/nextflow/bin/production_loop.py" "$@"
