#!/bin/bash
#
# cluster dialog within the mandelbrot project
# by R. Fruth, S. Kaiser & E. Kuhnt
#

MPI_PATH="/cvos/shared/apps/mpich2/ge/gcc/64/1.3.1/bin/mpirun"
BINARY="mandelbrot"

# checks if the binary is present
function checkBinary {
	# check if the binary exists
	if [ ! -f "$BINARY" ] ; then
		echo "ERROR: can not find mandelbrot binary! Make sure to run 'make' before using this dialog!"
		exit 1
	fi
}

# submit or output a job
function outputJob {

	echo "Do you want to submit your job (type: 0) or save it to a file (type: 1)?"
	read YESNO
	if [ "${YESNO}" == "0" ] ; then
		echo "Submitting your job..."
		echo -e "${1}" | qsub
	elif [ "${YESNO}" == "1" ] ; then
		echo "Please state the name ofthe file your job shall be saved to (default: job.pbs): "		
		read JOBFILE
		echo "Saving your job to file..."
		if [ ! "${JOBFILE}" == "" ] ; then
			echo -e "${1}" > ${JOBFILE}
		else
			echo -e "${1}" > "job.pbs"
		fi
	else
		echo "Oooops I did not understand you. But I assume you want to submit your job!"
		echo "submitting your job..."
		echo -e "${1}" | qsub
	fi

}

# prints the dialog header
function printDialogHeader {

	echo ""
	echo "mandelbrot job configurator"
	echo "(c) 2012 by R. Fruth, S. Kaiser & E. Kuhnt\n\n"
	echo -e "Welcome to the configuration of your mandelbrot job!\n\n"

}

# prints an introduction to the dialog
function printDialogIntro {

	echo "this dialog will help you performing the following configuration steps:"
	echo "1. configure your PBS environment (only the needed basics)"
	echo "2. configure your mandelbrot call"
	echo -e "3. submit your job or save the job file\n\n"

}

# main routine
function main {
	
	checkBinary
	
	printDialogHeader

	printDialogIntro

	echo "1. configure your PBS environment (only the needed basics)"
	# can not move them to a function as the catching
	# of the variable blocks the output
	echo "Type the name of your job (default: mb): "
	read PBS_NAME
	if [ "${PBS_NAME}" == "" ] ; then
		PBS_NAME="#PBS -N mb"
	else
		PBS_NAME="#PBS -N ${PBS_NAME}"
	fi

	echo "Do you want to get notified via eMail if your job exits or aborts? Type your mail address: "
	read PBS_MAIL
	if [ ! "${PBS_MAIL}" == "" ] ; then
		PBS_MAIL="#PBS -m ae\n#PBS -M ${PBS_MAIL}"
	fi
	
	echo "How many nodes do your want to reserve (default: 4)? "
	read PBS_NODES
	if [ "${PBS_NODES}" == "" ] ; then
		PBS_NODES="#PBS -l nodes=4"
	else
		PBS_NODES="#PBS -l nodes=${PBS_NODES}"
	fi

	echo "How many cores per node do you want to reserve (default: 1)? "
	read PBS_PPN
	if [ "${PBS_PPN}" == "" ] ; then
		PBS_PPN=":ppn=1"
	else
		PBS_PPN=":ppn=${PBS_PPN}"
	fi

	echo "How much walltime do you want to reserve (default: 00:00:05)?"
	read PBS_WALLTIME
	if [ "${PBS_WALLTIME}" == "" ] ; then
		PBS_WALLTIME="#PBS -l walltime=00:00:05"
	else
		PBS_WALLTIME="#PBS -l walltime=${PBS_WALLTIME}"
	fi

	echo -e "Please set your mpi path\n(default: ${MPI_PATH}): "
	read MPI_PATH
	if [ "${MPI_PATH}" == "" ] ; then
		MPI_PATH="/cvos/shared/apps/mpich2/ge/gcc/64/1.3.1/bin/mpirun"
	fi

	echo "How many slaves shall mpi create (default: 4)? "
	read MPI_NODES
	if [ "${MPI_NODES}" == "" ] ; then
		MPI_NODES="-np 4"
	else
		MPI_NODES="-np ${MPI_NODES}"
	fi

	echo -e "\n2. configure your mandelbrot call"
	MB_ARGS=""

	echo "Do you want to use the provided best picture configuration (type: yes)?"
	read YESNO
	if [ "${YESNO}" == "yes" ] ; then
		MB_ARGS="--bestpic"
	else
		echo "You typed '${YESNO}' - so i assume you want to use your own configuration..."
		echo "What shall be the name of your bitmap file (default: 'mandelbrot.bmp')? "
		read MB_FILENAME
		if [ ! "${MB_FILENAME}" == "" ] ; then
			MB_FILENAME="--filename ${MB_FILENAME}"
		fi

		echo "What shall be the width of your picture (default: 800)? "
		read MB_WIDTH
		if [ ! "${MB_WIDTH}" == "" ] ; then
			MB_WIDTH="--width ${MB_WIDTH}"
		fi

		echo "What shall be the height of your picture (default: 600)? "
		read MB_HEIGHT
		if [ ! "${MB_HEIGHT}" == "" ] ; then
			MB_HEIGHT="--height ${MB_HEIGHT}"
		fi

		echo "How many time shall be iterated over a point (default: 120)? "
		read MB_ITERATIONS
		if [ ! "${MB_ITERATIONS}" == "" ] ; then
			MB_ITERATIONS="-n ${MB_ITERATIONS}"
		fi

		echo "How many slices shall the picture be splitted to (default: 11)? "
		read MB_SLICES
		if [ ! "${MB_SLICES}" == "" ] ; then
			MB_SLICES="-s ${MB_SLICES}"
		fi

		echo "What shall be the red value of the mandelbrot set (default: 255)? "
		read MB_RED
		if [ ! "${MB_RED}" == "" ] ; then
			MB_RED="--red ${MB_RED}"
		fi

		echo "What shall be the green value of the mandelbrot set (default: 255)? "
		read MB_GREEN
		if [ ! "${MB_GREEN}" == "" ] ; then
			MB_GREEN="--green ${MB_GREEN}"
		fi

		echo "What shall be the blue value of the mandelbrot set (default: 0)? "
		read MB_BLUE
		if [ ! "${MB_BLUE}" == "" ] ; then
			MB_BLUE="--blue ${MB_BLUE}"
		fi

		echo "What shall be the background red value (default: 0)? "
		read MB_REDBG
		if [ ! "${MB_REDBG}" == "" ] ; then
			MB_REDBG="--redBg ${MB_REDBG}"
		fi

		echo "What shall be the background green value (default: 0)? "
		read MB_GREENBG
		if [ ! "${MB_GREENBG}" == "" ] ; then
			MB_GREENBG="--greenBg ${MB_GREENBG}"
		fi

		echo "What shall be the background blue value (default: 0)? "
		read MB_BLUEBG
		if [ ! "${MB_BLUEBG}" == "" ] ; then
			MB_BLUEBG="--blueBg ${MB_BLUEBG}"
		fi

		echo "What shall be the maximal real value (default: 1.0)? "
		read MB_REMAX
		if [ ! "${MB_REMAX}" == "" ] ; then
			MB_REMAX="--remax ${MB_REMAX}"
		fi

		echo "What shall be the minimal real value (default: -2.0)? "
		read MB_REMIN
		if [ ! "${MB_REMIN}" == "" ] ; then
			MB_REMIN="--remin ${MB_REMIN}"
		fi

		echo "What shall be the minimal imaginary value (default: -1.2)? "
		read MB_IMMIN
		if [ ! "${MB_IMMIN}" == "" ] ; then
			MB_IMMIN="--immin ${MB_IMMIN}"
		fi

	
		MB_ARGS="${MB_FILENAME} ${MB_WIDTH} ${MB_HEIGHT} ${MB_ITERATIONS} ${MB_SLICES}"
		MB_ARGS="${MB_ARGS} ${MB_RED} ${MB_GREEN} ${MB_BLUE} ${MB_REDBG} ${MB_GREENBG} ${MB_BLUEBG}"
		MB_ARGS="${MB_ARGS} ${MB_REMAX} ${MB_REMIN} ${MB_IMMIN}"
	fi

	PBS_CONFIG="${PBS_NAME}\n#PBS -j oe\n#PBS -r n\n${PBS_MAIL}\n${PBS_WALLTIME}"
	PBS_CONFIG="${PBS_CONFIG}\n${PBS_NODES}${PBS_PPN}\ncd \$PBS_O_WORKDIR"
	PBS_CONFIG="${PBS_CONFIG}\n${MPI_PATH} ${MPI_NODES} -machinefile \$PBS_NODEFILE `pwd`/${BINARY}"
	JOB="#!/bin/bash\n${PBS_CONFIG} ${MB_ARGS}"
	
	echo -e "\n3. submit your job or save the job file"
	outputJob "$JOB"
	
}

main "$@"
