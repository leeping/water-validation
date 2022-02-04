#!/bin/bash

# This wrapper script is for running GROMACS jobs on clusters.

# Load my environment variables. :)
. /etc/profile
. ~/.bashrc

# Disable GROMACS backup files
export GMX_MAXBACKUP=-1

export PATH=$HOME/opt/gromacs/bin:$PATH
export LD_LIBRARY_PATH=$HOME/opt/gromacs/lib:$LD_LIBRARY_PATH
if [[ $HOSTNAME =~ "leeping" ]] ; then
    . /opt/intel/Compiler/11.1/072/bin/iccvars.sh intel64
elif [[ $HOSTNAME =~ "sh" ]] ; then
    . /opt/intel/composerxe/bin/compilervars.sh intel64
elif [[ $HOSTNAME =~ "cn" ]] ; then
    . $HOME/opt/intel/bin/compilervars.sh intel64
elif [[ $HOSTNAME =~ "vsp" ]] ; then
    . $HOME/opt/intel/bin/compilervars.sh intel64
elif [[ $HOSTNAME =~ "biox" ]] ; then
    . $HOME/opt/intel/bin/compilervars.sh intel64
fi

echo "#=======================#"
echo "# ENVIRONMENT VARIABLES #"
echo "#=======================#"
echo

env

echo
echo "#=======================#"
echo "# STARTING CALCULATION! #"
echo "#=======================#"
echo

nt=1
if [ ! -z $OMP_NUM_THREADS ] ; then nt=$OMP_NUM_THREADS ; fi

# Gromacs looks for this environment variable and behaves weirdly
unset OMP_NUM_THREADS

grompp -f eq.mdp -o eq.tpr -c min.gro -maxwarn 2 2>&1
mdrun -v -nt $nt -deffnm eq 2>&1
grompp -c eq.gro -f md.mdp -o md.tpr -maxwarn 2 2>&1
mdrun -v -nt $nt -deffnm md 2>&1

for i in "eq.gro" "eq.tpr" "eq.xtc" "eq.trr" "eq.edr" "md.gro" "md.tpr" "md.xtc" "md.trr" "md.edr" ; do
    if [ -f $i ] ; then
        tar rvf result.tar $i
    fi
done

bzip2 result.tar
