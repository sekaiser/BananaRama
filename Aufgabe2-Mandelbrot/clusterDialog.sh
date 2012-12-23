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

	echo "Type the name of your job (default: mb): "
	read PBS_NAME
	if [ "$PBS_NAME" == "" ] ; then
		PBS_NAME="#PBS -N mb"
	else
		PBS_NAME="#PBS -N $PBS_NAME"
	fi

	echo "Do you want to get notified via eMail if your job exits or aborts? Type your mail address: "
	read PBS_MAIL
	if [ ! "$PBS_MAIL" == "" ] ; then
		PBS_MAIL="#PBS -m ae\n#PBS -M $PBS_MAIL"
	fi
	
	PBS_MAILADDR=""		# user mail address
	PBS_WALLTIME=""		# wall time on the cluster
	PBS_NODES=""		# number of nodes
	PBS_PPN=""		# cores per node

	MPI_PATH="/cvos/shared/apps/mpich2/ge/gcc/64/1.3.1/bin/mpirun" # mpi path
	MPI_NODES="-np 4"

	echo -e "PBSmail: $PBS_MAIL"

}

main "$@"
