#!/bin/bash
#SBATCH --partition=pre
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=gHg_gromacs_test

module load gcc cmake openmpi
export LD_LIBRARY_PATH=$OPENMPI_ROOT/lib:$LD_LIBRARY_PATH

ulimit -s unlimited
export OMP_NUM_THREADS=16
export OMP_STACKSIZE=256M
export PYSCF_MAX_MEMORY=64000

fnms='gHg'

# Check if SLURM_JOB_ID is set
if [ -z "$SLURM_JOB_ID" ]; then
    echo "Error: SLURM_JOB_ID is not set. Are you running this script within a SLURM job?" >&2
    exit 1
fi

# Base directory for temporary folders
BASE_DIR="/local/$USER/$SLURM_JOB_ID"

# Ensure the base directory exists
mkdir -p "$BASE_DIR"

# Function to run GROMACS with a specific PYSCF_TMPDIR
run_gromacs() {
    local tmp_dir=$1
    export PYSCF_TMPDIR="$tmp_dir"
    fnm=$2
    outfile=${fnm}.qmmmout
    echo "Starting" > $outfile
    date >> $outfile
        echo 'running on node: '$SLURM_JOB_NODENAME > $outfile
    echo   'Number of nodes allocated                            : '  $SLURM_JOB_NUM_NODES       >> $outfile
    echo   'Total number of tasks/processes requested (--ntasks) : '  $SLURM_NTASKS              >> $outfile
    echo   'Tasks per node                                       : '  $SLURM_NTASKS_PER_NODE     >> $outfile
    echo   'Total CPUs on the allocated node                     : '  $SLURM_CPUS_ON_NODE        >> $outfile
    echo   'allocated per task (--cpus-per-task)                 : '  $SLURM_CPUS_PER_TASK	CPUs >> $outfile
    echo   'CPUs per node allocated to the job                   : '  $SLURM_JOB_CPUS_PER_NODE   >> $outfile
    echo   'List of nodes assigned                               : '  $SLURM_NODELIST            >> $outfile
    echo   'Local task ID on the node                            : '  $SLURM_LOCALID             >> $outfile
    echo   'Global MPI rank (0-based task ID)                    : '  $SLURM_PROCID              >> $outfile
    echo   'Synonym for SLURM_NTASKS                             : '  $SLURM_NPROCS              >> $outfile
    sed -n '10,39p;40q' pyscfdriver.py  >> $outfile
    srun --exclusive -N 1 -n 1 -c 16 bash -c "mkdir -p $tmp_dir && ~/github/gromacs/build/bin/gmx_d mdrun -deffnm $fnm " >> $outfile 2>&1

}

for i in $fnms; do
    # if [ ! -f "$i.gro" ]; then
        tmp_dir="$BASE_DIR/$i"
        # mkdir -p "$tmp_dir"  # Ensure each tmp directory exists -> This is moved to srun command
        (run_gromacs "$tmp_dir" "$i") &
    # fi
done

# Wait for all background processes to finish
wait

rm -rf "$BASE_DIR"
