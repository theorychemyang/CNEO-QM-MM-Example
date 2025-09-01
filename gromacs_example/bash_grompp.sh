tnm=gHg
logfile=${tnm}.gromppout
rm $logfile
cnm=../conf
fnm=./run_params
nnm=../index
pnm=../topol

rm $tnm.*
date > $logfile

$HOME/github/gromacs/build/bin/gmx_d grompp -f $fnm.mdp -c $cnm.gro -p $pnm.top -n $nnm.ndx -o ${tnm}.tpr -maxwarn 10 >> $logfile 2>&1
cat $fnm.mdp >> $logfile

rm \#$tnm*
