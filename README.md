#External Python User Interface for MCCCS Towhee

###Notes
The two Python scripts in source have been written in Python3.5 and are intended as a User Interface (UI) for 
isobaric-isothermal, canonical, and Gibbs ensemble Monte Carlo simulations with Monte Carlo for Complex 
Chemical Systems (MCCCS) Towhee version 7.1.0.  

The software has been compiled with Python3.5 using the PyCharm Community Edition, and with Python3.4 using
the Ubuntu command line. No warnings or weak warnings exist in the current version of the code.  

###Details
The UI has been designed as an organziational, time saving, and data analysis tool for performing molecular 
simulations with Towhee on both personal computers, and clusters. The UI is composed of two scripts 
- towhee_initialize.py and towhee_collect.py - and has the following optional features: 

1) towhee_intialize.py
  - creates pressure (P), temperature(T), and set specific directories.
  - automatically copy 'towhee_input' files at specified T and P to 
    each final directory with different random number seeds
  - automatically start simulations for equilibration and production
    stages
  - automatically copy intial system configuration and molecule template
    files to directories
  - specify numerous input specifications from http://towhee.sourceforge.net/
    input/towhee_input_v7_0_x.html within script using a generic 
    'towhee_input' file or specify all input in the 'towhee_input' file 
  - automatically generate necessary 'towhee_parallel' and MPI executable 
    files within T directories for simulations on a cluster
    
2) towhee_collect.py
  - automatically collect various thermodynamic property results from output
    files from a single T and P for any number of sets, compute averages, 
    standard deviations, and percent standard deviations, and print all data 
    to an output file

###Contributors
*Robert Kelly
