To compile de program: gfortran -o dynamics dynamics.f90 Then, a file called 'dynamics' will be created. To run the program: ./dynamics 
In this moment it will ask you: 'Write the name of your input file' You write: ar.xyz (or the name of the file that contains the coordinates in the required format) 
Then the program will ask you: 'Each how many iterations you want to write info:' You need to introduce between how many steps you want to have info about 
Then the program will run and a file called trajectory.xyz is created where the info about the dynamics is stored.
