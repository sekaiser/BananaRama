#!/bin/bash
#
# cluster dialog within the mandelbrot project
# by R. Fruth, S. Kaiser & E. Kuhnt
#

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

}

main "$@"
