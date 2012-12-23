#!/bin/bash
#
# cluster dialog within the mandelbrot project
# by R. Fruth, S. Kaiser & E. Kuhnt
#

MB_WIDTH=""
MB_HEIGHT=""
MB_ITERATIONS=""
MB_SLICES=""
MB_BESTPIC=""
MB_FILENAME=""
MB_RED=""
MB_GREEN=""
MB_BLUE=""
MB_REDBG=""
MB_GREENBG=""
MB_BLUEBG=""
MB_REMAX=""
MB_REMIN=""
MB_IMMIN=""

# checks if the binary is present
function checkBinary {
	# check if the binary exists
	if [ ! -f mandelbrot ] ; then
		echo "ERROR: can not find mandelbrot binary! Make sure to run 'make' before using this dialog!"
		exit 1
	fi
}

# prints the dialog header
function printDialogHeader {

	echo ""
	echo "mandelbrot job configurator"
	echo "(c) 2012 by R. Fruth, S. Kaiser & E. Kuhnt"
	echo ""
	echo "Welcome to the configuration of your mandelbrot job!"
	echo ""
	echo ""

}

# prints an introduction to the dialog
function printDialogIntro {

	echo "this dialog will help you performing the following configuration steps:"
	echo "1. configure your PBS environment (only the needed basics)"
	echo "2. configure your mandelbrot call"
	echo "3. submit your job or save the job file"

}

# main routine
function main {
	
	checkBinary
	
	printDialogHeader

	printDialogIntro

	# can not move them to a function as the catching
	# of the variable blocks the output
	echo "Type the name of your job (default: mb): "
	read PBS_NAME
	if [ "$PBS_NAME" == "" ] ; then
		PBS_NAME="#PBS -N mb\n"
	else
		PBS_NAME="#PBS -N $PBS_NAME'n"
	fi

	echo "Do you want to get notified via eMail if your job exits or aborts? Type your mail address: "
	read PBS_MAIL
	if [ ! "$PBS_MAIL" == "" ] ; then
		PBS_MAIL="#PBS -m ae\n#PBS -M $PBS_MAIL\n"
	fi
	
	echo "How many nodes do your want to reserve (default: 4)? "
	read PBS_NODES
	if [ "$PBS_NODES" == "" ] ; then
		PBS_NODES="#PBS -l nodes=4"
	else
		PBS_NODES="#PBS -l nodes=$PBS_NODES"
	fi

	echo "How many cores per node do you want to reserve (default: 1)? "
	read PBS_PPN
	if [ "$PBS_PPN" == "" ] ; then
		PBS_PPN=":ppn=1\n"
	else
		PBS_PPN=":ppn=$PBS_PPN\n"
	fi

	echo "How much walltime do you want to reserve (default: 00:00:05)?"
	read PBS_WALLTIME
	if [ "$PBS_WALLTIME" == "" ] ; then
		PBS_WALLTIME="#PBS -l walltime=00:00:05"
	else
		PBS_WALLTIME="#PBS -l walltime=$PBS_WALLTIME\n"
	fi



	
	MPI_PATH="/cvos/shared/apps/mpich2/ge/gcc/64/1.3.1/bin/mpirun"


	MPI_NODES="-np 4"

	echo -e ""

}

main "$@"
