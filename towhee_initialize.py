# Python 3.4.0/3.5.0 -- Towhee v7.1.0

# This program requires three command and function libraries: the 'os' for operating system commands, 'shutil'
# for copy commands, and 'random' to generate random numbers.

import os
import shutil
import random

# Start towhee_initialize

# The function to be performed; depends on if an initial system configuration is available or if a
# configuration must be generated.

function = "equilibration"
# function = "production"

# The platform that the simulations will be performed on. This string allows for logical creation of files
# that are required on a cluster but are not required on a personal computer.

platform = "pc"
# platform = "cluster"

# These logical variables specify whether or not to auto create directories, auto copy input files, auto start
# simulations, and/or auto copy template and initial configuration files into set directories. Note that there is no
# auto start feature for cluster processing.

createDir = True
autoCopy = True
autoStart = True
template = False
config = False

# The species directory 'speciesDir' specifies the name of the system or mixture system that the main directory will
# be created as. The 'ensemble' is used to define the ensemble in 'towhee_input'. Set, temperature (K), and pressure
# (bar) lists specify the number of sets at each temperature, the temperatures at
# each pressure, and each pressure. 'numSims' is the number of simulations that will be performed.

speciesDir = "methane"
ensemble = "npt"
initDir = os.getcwd()
setList = ["set1", "set2", "set3", "set4"]
tempList = ["300.0"]
pressList = ["50.0"]
numSims = int(len(setList)) * int(len(tempList)) * int(len(pressList))

# 'towheeInput' specifies the name of the input file that will be copied into each set directory

towheeInput = "towhee_input"

# 'user' and 'eMail' allow the user to specify the user name on the cluster and the users email for use in
# creation of necessary files. Note that 'user' is also the part of the path to force field files
# and to the working directory that 'towhee_initialize.py' will be executed from on clusters - this assumes
# that for clusters the program is run from the directory '/home/user/'. 'eMail' is the email for situations
# where an email client is setup.

user = "robert"
eMail = "user@email.com"

# 'nMolec' is the number of molecules of each species in the system - read by towhee as a one-dimensional
# array with an element for each molecule type; 'nStep' is the number of simulation steps; 'nMolecType'
# is the number of different types of molecules in the system

nMolec = str("100")
nStep = 10000
nMolecType = 1

# 'numBoxes' is the number of simulation boxes; 'printFreq' is the step interval before the energy is printed;
# 'blockSize' defines the number of blocks used for intra-simulation averaging; 'movieFreq' is the frequency
# for storing system configurations in the 'towhee_movie' file; 'backupFreq' is the frequency that towhee backup
# files are created in each set directory; 'restartFreq' is the frequency that towhee restart files are created in
# each set directory; 'pbdFreq' is the frequency that .pbd files are created; 'pressFreq' is the frequency for
# calculating the virial pressure; 'chemPot' is the number of trial insertions for each box. Allows for calculation
# of the chemical potential in non-grandcanonical simulations.

numBoxes = 1
printFreq = int(nStep / 10)
blockSize = int(nStep / 10)
movieFreq = 0
backupFreq = int(nStep / 4)
restartFreq = int(nStep / 4)
pbdFreq = int(nStep / 4)
pressFreq = 0
chemPot = str("0")

# 'numFF' is the number of force fields; 'forceField' is a the potential energy force field; 'pathFF' is a
# concatenated string (includes user) for the path to the force field file. Note that the concatenated string
# for the force field path assumes all towhee force fields are found in '/home/user/ForceFields'.

numFF = 1
forceField = "TraPPE-UA"
pathFF = "/home/" + user + "/ForceFields/" + "towhee_ff_" + forceField

# 'potential' is the potential function used to calculate the system energy; 'mixRule' is the mixing rule used in
# the potential energy function; 'rMin' is the minimum intermolecular distance between atoms in the simulation;
# 'rCut' is the radial cutoff distance for the potential function; 'rCutin' is the inner nonbonded cutoff for
# configurational bias moves; 'initStyle' specifies the initial box style; 'initLattice' is the type of lattice;
# 'initMol' specifies the initial number of molecules of each type; 'initXYZ' specifies the number of
# three-dimensional grid points for the initial configuration; 'dimX', 'dimY', and 'dimZ' are the initial box
# dimensions; 'trmaxFreq' is the frequency that the max translational displacement is updated; 'volmaxFreq' is
# the frequency that the max volume displacement is updated; 'electrostatics' is a logical to specify
# electrostatics, and the input file details for Coulombic interactions; 'chargeAssign' is the method of charge
# assignment to be performed by towhee (i.e., bond-increment, manual, etc.).

potential = "Lennard-Jones"
mixRule = "Lorentz-Berthelot"
rMin = 1.0
rCut = 7.0
rCutin = 7.0
initStyle = str("full cbmc")
initLattice = str("simple cubic")
initMol = str("100")
initXYZ = "5 5 4"
dimX = "20.0D0 0.00D0 0.00D0"
dimY = "0.00D0 20.0D0 0.00D0"
dimZ = "0.00D0 0.00D0 20.0D0"
trmaxFreq = "trmaxFreq"
volmaxFreq = "volmaxFreq"
electroStatics = "none"
chargeAssign = "none"

# If the force field requires electrostatics, standard towhee treatment of Coulombic interactions are written to the
# input file. 'charge' is the character string for the electrostatic specifications.

charge = "charge"
if electroStatics == "coulomb":
    charge = "'" + "coulomb" + "'" + "\n" \
             "coulombstyle" + "\n" \
             "'" + "ewald_fixed_kmax" + "'" + "\n" \
             "kalp" + "\n" \
             "5.6d0" + "\n" \
             "kmax" + "\n" \
             "5" + "\n" \
             "dielect" + "\n" \
             "1.0d0"

elif electroStatics == "none":
    charge = "'" + "none" + "'"

# Depending on the function ("equilibration" or "production") and whether an initial configuration will be used,
# trmaxFreq, volmaxFreq, and lInit are treated differently.  This is treated with four 'if' statements below.

lInit = "lInit"
if function == "equilibration" and template is False and config is False:
    trmaxFreq = int(nStep / 1000)
    volmaxFreq = int(nStep / 1000)
    lInit = ".true."

elif function == "equilibration" and template is True and config is False:
    trmaxFreq = int(nStep / 1000)
    volmaxFreq = int(nStep / 1000)
    lInit = ".true."
    initStyle = str("template")

elif function == "equilibration" and template is False and config is True:
    trmaxFreq = int(nStep / 1000)
    volmaxFreq = int(nStep / 1000)
    lInit = ".false."

elif function == "production":
    trmaxFreq = int(nStep / 10)
    volmaxFreq = int(nStep / 10)
    lInit = ".false."

# The generic or 'base' input file from the initial directory is opened and all input data is read as 'allInput'

file = open(towheeInput)
allInput = file.read()
file.close()

# A triple 'for' loop in pressure, temperature, and sets is performed to concatenate strings to the final directory
# based on the pressures in pressList, temperatures in tempList, and sets in setList. After the pressure directory
# is created, the pressure is converted from bar (my preferred units) to kPa (pressure units in towhee).

for P in pressList:
    pressDir = speciesDir + "/" + P
    P = float(P) * float(100.0)
    P = str(P)

    for T in tempList:
        tempDir = pressDir + "/" + T

        for Set in setList:
            finalDir = tempDir + "/" + Set + "/"

            # If auto create directories is specified, all directories are created using the os.makedirs() function and
            # the concatenated final directory string(s).

            if createDir is True:
                os.makedirs(finalDir)

            # If a towhee_template file is specified to be copied into each final directory, the shutil copy function
            # is used to copy "towhee_template" to each final directory. Note this file must be in the same directory
            # as 'towhee_initialize.py'.

            if template is True and config is False and function == "equilibration":
                shutil.copy("towhee_template", finalDir)

            # If a towhee_template is not desired as an initial configuration for each molecule, and an initial system
            # configuration is specified, the shutil copy function is used to copy "towhee_initial" to each final
            # directory. Again, 'towhee_initial' must be in the same working directory as 'towhee_initialize.py'

            if config is True and template is False:
                shutil.copy("towhee_initial", finalDir)

            # From the base input file that has been read to allInput, specific character strings are replaced with
            # the numerical values and strings described/specified earlier.

            specificInput = allInput.replace("insert_random", str(random.randint(0, 100000)))
            specificInput = specificInput.replace("ensemble_input", str("'" + ensemble + "'"))
            specificInput = specificInput.replace("temperature_input", str(T + "D0"))
            specificInput = specificInput.replace("pressure_input", str(P + "D0"))
            specificInput = specificInput.replace("nmolty_input", str(nMolecType))
            specificInput = specificInput.replace("nmolec_input", str(nMolec))
            specificInput = specificInput.replace("numboxes_input", str(numBoxes))
            specificInput = specificInput.replace("nstep_input", str(nStep))
            specificInput = specificInput.replace("print_input", str(printFreq))
            specificInput = specificInput.replace("blocksize_input", str(blockSize))
            specificInput = specificInput.replace("moviefreq_input", str(movieFreq))
            specificInput = specificInput.replace("backup_input", str(backupFreq))
            specificInput = specificInput.replace("restart_input", str(restartFreq))
            specificInput = specificInput.replace("pbd_input", str(pbdFreq))
            specificInput = specificInput.replace("pressure_freq", str(pressFreq))
            specificInput = specificInput.replace("trmax_freq", str(trmaxFreq))
            specificInput = specificInput.replace("volmax_freq", str(volmaxFreq))
            specificInput = specificInput.replace("chempot_per_step", str(chemPot))
            specificInput = specificInput.replace("numff_input", str(numFF))
            specificInput = specificInput.replace("ff_input", str(pathFF))
            specificInput = specificInput.replace("potential_input", str("'" + str(potential) + "'"))
            specificInput = specificInput.replace("mixrule_input", str("'" + str(mixRule) + "'"))
            specificInput = specificInput.replace("rmin_input", str(str(rMin) + "D0"))
            specificInput = specificInput.replace("rcut_input", str(str(rCut) + "D0"))
            specificInput = specificInput.replace("rcutin_input", str(str(rCutin) + "D0"))
            specificInput = specificInput.replace("electro_input", str(charge))
            specificInput = specificInput.replace("linit_input", str(lInit))
            specificInput = specificInput.replace("initstyle_input", str("'" + str(initStyle) + "'"))
            specificInput = specificInput.replace("initlattice_input", str("'" + str(initLattice) + "'"))
            specificInput = specificInput.replace("initmol_input", str(initMol))
            specificInput = specificInput.replace("initxyz_input", str(initXYZ))
            specificInput = specificInput.replace("dimx_input", str(dimX))
            specificInput = specificInput.replace("dimy_input", str(dimY))
            specificInput = specificInput.replace("dimz_input", str(dimZ))
            specificInput = specificInput.replace("ff2_input", str("'" + str(forceField) + "'"))
            specificInput = specificInput.replace("charge_input", str("'" + str(chargeAssign) + "'"))

            # Specific functions are carried out for equilibration simulations

            if function == "equilibration":

                # If autoCopy is specified, 'specificInput' from above is copied into 'towhee_input' files in each final
                # directory.

                if autoCopy is True:
                    file = open(finalDir + towheeInput, 'w')
                    file.write(specificInput)
                    file.close()

                # If 'autoStart' is specified, the platform is a personal computer, and the number of simulations is
                # not more than 12 (assuming a quadcore pc - use number of processors * 3 to avoid autokilling
                # simulations), the simulations with be started automatically using the os.system() command line
                # function and the equilibration output directed to an output file.

                if numSims <= int(12) and autoStart is True and platform == "pc":
                    os.chdir(finalDir)
                    os.system("/usr/bin/nohup towhee > towhee_output1.txt &")
                    os.chdir(initDir)

            # If the function is production and config is false, the final configuration is copied to the initial
            # configuration, and other functions are carried out.

            if function == "production" and config is False:
                shutil.copy2(finalDir + "towhee_final", finalDir + "towhee_initial")

                # If autocopy is specified, 'specificInput' from above (and for production) is copied into
                # 'towhee_input' files in each final directory.

                if autoCopy is True:
                    file = open(finalDir + towheeInput, 'w')
                    file.write(specificInput)
                    file.close()

                # If 'autoStart' is specified, the platform is a personal computer, and the number of simulations is
                # not more than 12, the simulations with be started automatically using the os.system() command line
                # function and the production output directed to an output file.

                if numSims <= int(12) and autoStart is True and platform == "pc":
                    os.chdir(finalDir)
                    os.system("/usr/bin/nohup towhee > towhee_output2.txt &")
                    os.chdir(initDir)

# This portion of the code is specific for running simulations on a cluster, and creates additional files for the
# process manager and for 'towhee_parallel'.

# Depending on the type of simulation, different output files names are specified. 'towheeOutput' refers to the
# file towhee will print output data to. 'towheeReport' refers to the output file that the cluster process manager will
# print process information from the jobs submitted to a particular node.

towheeOutput = "towheeOutput"
if function == "equilibration":
    towheeOutput = "towhee_output1.txt"
    towheeReport = "towhee_report1"

elif function == "production":
    towheeOutput = "towhee_output2.txt"
    towheeReport = "towhee_report2"

# For cluster processing, parallel input is specified to be copied to a 'towhee_parallel' file in each temperature
# directory. 'numJobs' is the number of simulations to be processed at each temperature step.

if platform == "cluster":
    numJobs = str(int(len(setList)))
    parallelInput = "#number of jobs \n" + \
                    numJobs + "\n" \
                    "#stdout filename \n" + \
                    towheeOutput + "\n" \
                    "#working directories \n" \

    # A triple looped is performed again to specify processing information in the 'towhee_parallel' and 'base.batch'
    # files for cluster processing.

    for P in pressList:
        pressDir = speciesDir + "/" + P

        for T in tempList:
            tempDir = pressDir + "/" + T

            # The directory is changed to the temperature directory using the os.chdir() function, a
            # 'towhee_parallel' file is created, and parallelInput is copied into the file.

            os.chdir(initDir + "/" + tempDir)
            file = open("towhee_parallel", "w")
            file.write(parallelInput)
            file.close()

            # Information for the 'base.batch' file is specified as 'baseInput' using information specified earlier
            # and from the the double loop above.

            baseInput = "#!/bin/csh \n" \
                        "#$ -N " + speciesDir + "\n" \
                        "#$ -pe mpi " + numJobs + "\n" \
                        "#$ -e /home/" + user + "/" + tempDir + "/" + towheeReport + "\n" \
                        "#$ -o /home/" + user + "/" + tempDir + "/" + towheeReport + "\n" \
                        "#$ -cwd \n" \
                        "#$ -j y \n" \
                        "#$ -m be -M " + eMail + "\n" \
                        "cd /home/" + user + "/" + tempDir + "/" + "\n" \
                        "mpirun -np $NSLOTS -machinefile $TMPDIR/machines /home/" + user + "/bin/towhee" "\n" \

            # A 'base.batch' file is created, and 'baseInput' is copied into the file.

            file = open("base.batch", "w")
            file.write(baseInput)
            file.close()

            # The last loop of the triple loop is used to specify more information for the 'towhee_parallel' file. The
            # file is opened and the working directories for each set is printed to the file below the existing text.

            for Set in setList:
                finalDir = tempDir + "/" + Set + "/"
                workingDirs = "/home/" + user + "/" + finalDir

                os.chdir(initDir + "/" + tempDir)
                file = open("towhee_parallel", "a")
                file.write(workingDirs + "\n")
                file.close()

# End towhee_initialize
