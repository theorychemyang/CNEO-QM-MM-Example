# source /software/groups/chemistry_yang_group/cp2k/tools/toolchain/install/setup
# for i in 1.40 ; do
# nm=$i
tnm=gHg
logfile=${tnm}.gromppout
rm $logfile
cnm=./conf
# cnm=./r_2.55/r_2.55
# fnm=./f_qmmm_em_pull
fnm=./run_params_no_constraint
nnm=./index
pnm=./topol

rm $tnm.*
date > $logfile
# srun --mpi=pmix -N1 -n1 -p int --pty
# if [ ! -f $cnm.gro ]; then
#     echo "Error: $cnm.gro not found"
#     $HOME/github/gromacs/build/bin/gmx_d trjconv -f ./${tnm}/${tnm}.trr -s ./${tnm}/${tnm}.tpr -o $cnm.gro -dump 1000
# fi
$HOME/github/gromacs/build/bin/gmx_d grompp -f $fnm.mdp -c $cnm.gro -p $pnm.top -n $nnm.ndx -o ${tnm}.tpr -maxwarn 10 >> $logfile 2>&1

cat $fnm.mdp >> $logfile
rm \#$tnm*
# done
