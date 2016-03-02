# tasks_priority_matrix
Matlab scripts for testing a novel redundancy resolution method based on the tasks priority matrix

With the this code it possible to obtain the plots ans numerical results presented in the paper 
"Thoward the complete task state control" submitted to IROS 2016

The main mathod is reported in the function tasksPriorityMatrix.m.

To obtain the data in Table 1 please run findBigErr3T.m
Figure 1 is obtained with the script plotTvarNtasks.m and Figure 2 with plotTvarNdof.m

For the Romeo simulation:
run simulationRomeo_std for using the standard method or simulatingRomeo_tpm for using the presented method
then run videoSimulationRomeo.m to see the results
The results of the simulation are stored in workspace variable that can be used for statistics, please use them at your will.

