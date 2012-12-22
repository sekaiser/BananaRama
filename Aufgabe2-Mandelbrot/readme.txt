23.12.2012

content of this file
	- compiling the program
	- running the program (includes a local dialog)
	- running the cluster dialog
	- git

compiling the program
	Before compiling please make sure that mpich2, make, etc. are installed on
	your machine.

	To compile the program you just need to type 'make' within the project directory.
	The provided Makefile will place a mandelbrot binary in the project directory.

running the program
	After you compiled the program you may run it with the help parameter in order
	to see all available configuration parameters. The usage information itself is
	entirely self-explanatory and shouldn't leave you with open questions.

	You may print the usage information with the command:
		mpirun ./mandelbrot --help

running the program on the cluster
	Before running the cluser dialog make sure to compile the program!

	To run the program on the cluster you will either need a pbs file or use the
	clusterDialog.sh bash script, which will help you configure and run the 
	program.
	
	Since we provided you with a pretty comfy dialog, we did not give you a PBS file
	anyway the clusterDialog.sh bash script can generate one for you.

	You may start the clusterDialog.sh script with the command:
		bash clusterDialog.sh

git
	our git repo can be cloned from https://github.com/sekaiser/BananaRama.git

