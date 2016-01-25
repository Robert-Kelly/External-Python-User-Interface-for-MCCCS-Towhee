# Python 3.4.0/3.5.0 -- Towhee v7.1.0

# This program requires two command and function libraries: the 'os' for operating system commands, and 'statistics'
# for computing towhee output statistics.

import os
import statistics

# Start towhee_collect

# The species directory specifies the name of the initial directory. Set, temperature (K), and pressure (bar) lists
# specify the sets at each temperature, the temperatures at each pressure, and each pressure. This information
# is used to concatenate strings in the same way as towhee_initialize.

# Note that this script can only gather towhee output data from a single temperature and pressure point for any
# number of sets.

speciesDir = "methane"
setList = ["set1", "set2", "set3", "set4"]
tempList = ["300.0"]
pressList = ["50.0"]

# 'ensemble' is the simulation type specifying the three constant state variables.  This variable functions as 'gibbs'
# for two simulation boxes, and arbitrarily as "npt" or "nvt" for single box simualations.

ensemble = "npt"

# 'function' specifies the simulation stage that has already been processed, and allows for assignment of the towhee
# output file to be read. Note that the output file names automatically corresponds to the names from
# towhee_initialize.

function = "equilibration"
# function = "production"

# The arrays defined below are used to append numerical data from each set output file. 'nMolTyList' are the number of
# molecule types from each simulation; 'nMolecTypList' are the FORTRAN two-dimensional array values for the number
# of molecules of each type;  'temperatureList' and 'pressureList' are the temperatures and pressures in each set;
# nMolecList are the total number of molecules from each set;  'molecMassList1' and 'molecMassList2' are the
# molecular mass values for a binary mixure simualtion;  'internalEnergyList' are the internal energies of departure
# from each set; 'molarVolumeList' are the molar volumes from each set.

nMolTyList = []
nMolecTypList = []
temperatureList = []
pressureList = []
nMolecList = []
molecMassList1 = []
molecMassList2 = []
internalEnergyList = []
molarVolumeList = []

# 'molecNumBox1' and 'molecNumBox2' are the number of molecules in box one and box two;  'molarVolumeBox1' and
# 'molarVolumeBox2' are the molar volumes in box one and two; 'virialPressBox1' and 'virialPressBox2' are the virial
# pressures in box one and two; 'thermoPressBox1' and 'thermoPressBox2' are the thermodynamic pressures in box
# one and two; 'chemPotBox1' and 'chemPotBox2' are the Gibbs chemical potentials in box one and box two; 'deltaHList'
# is the enthalpy of vaporization list.

molecNumBox1 = []
molecNumBox2 = []
molarVolumeBox1 = []
molarVolumeBox2 = []
virialPressBox1 = []
virialPressBox2 = []
thermoPressBox1 = []
thermoPressBox2 = []
chemPotBox1 = []
chemPotBox2 = []
deltaHList = []

# 'if' statements are used to define the file to be read depending on if the simulation stage is equilibration or
# production.

if function == "equilibration":
    towheeOutput = "towhee_output1.txt"

elif function == "production":
    towheeOutput = "towhee_output2.txt"

# A triple 'for' loop is performed to concatenate the set directory string(s) in the same way as towhee_initialize.
# The path to the current directory is defined for later use.

currentPath = str(os.getcwd() + "/")
for P in pressList:
    pressDir = currentPath + speciesDir + "/" + P

    for T in tempList:
        tempDir = pressDir + "/" + T

        for Set in setList:
            finalDir = tempDir + "/" + Set
            os.chdir(finalDir)

            # The towhee output file in each set is opened and all output is read to 'allOutput'

            file = open(towheeOutput)
            allOutput = file.readlines()
            file.close()

            # A loop is performed through all lines of text in 'allOutput' and if a specified string is found in the
            # line, the string is replace with nothing ("") leaving only the numerical value that corresponds to
            # the string.

            for line in allOutput:

                if ensemble != "gibbs":

                    # Specified strings are replaced for each array defined earlier. The left over numerical value from
                    # the string is appended to the defined array.

                    if line.find("nmolty:") != -1:
                        nMolTy = line.replace("nmolty:", "")
                        nMolTy = nMolTy.strip()
                        nMolTyList.append(int(nMolTy))

                    if line.find("nmolectyp:") != -1:
                        nMolecTyp = line.replace("nmolectyp:", "")
                        nMolecTyp = nMolecTyp.strip()
                        nMolecTypList.append(str(nMolecTyp))

                    if line.find("Temperature [K]:") != -1:
                        temperature = line.replace("Temperature [K]:", "")
                        temperature = temperature.strip()
                        temperatureList.append(float(temperature))

                    # Pressure is converted from kPa back to bar.

                    if line.find("External pressure [kPa]:") != -1:
                        pressure = line.replace("External pressure [kPa]:", "")
                        pressure = pressure.strip()
                        pressure = float(pressure) / 100
                        pressureList.append(float(pressure))

                    if line.find("Number of molecules:") != -1:
                        nMolec = line.replace("Number of molecules:", "")
                        nMolec = nMolec.strip()
                        nMolecList.append(int(nMolec))

                    if line.find("Molecular mass for molecule type     1 is") != -1:
                        molecMass1 = line.replace("Molecular mass for molecule type     1 is", "")
                        molecMass1 = molecMass1.replace("g/mol", "")
                        molecMass1 = molecMass1.strip()
                        molecMassList1.append(float(molecMass1))

                    if line.find("Molecular mass for molecule type     2 is") != -1:
                        molecMass2 = line.replace("Molecular mass for molecule type     2 is", "")
                        molecMass2 = molecMass2.replace("g/mol", "")
                        molecMass2 = molecMass2.strip()
                        molecMassList2.append(float(molecMass2))

                    if line.find("U                    kJ/mol") != -1:
                        internalEnergy = line.replace("U                    kJ/mol", "")
                        internalEnergy = internalEnergy.strip()
                        internalEnergyList.append(float(internalEnergy))

                    if line.find("Molar Volume         ml/mol") != -1:
                        molarVolume = line.replace("Molar Volume         ml/mol", "")
                        molarVolume = molarVolume.strip()
                        molarVolumeList.append(float(molarVolume))

                if ensemble == "gibbs":

                    if line.find("nmolty:") != -1:
                        nMolTy = line.replace("nmolty:", "")
                        nMolTy = nMolTy.strip()
                        nMolTyList.append(int(nMolTy))

                    if line.find("Temperature [K]:") != -1:
                        temperature = line.replace("Temperature [K]:", "")
                        temperature = temperature.strip()
                        temperatureList.append(float(temperature))

                    if line.find("Number of molecules:") != -1:
                        nMolec = line.replace("Number of molecules:", "")
                        nMolec = nMolec.strip()
                        nMolecList.append(int(nMolec))

                    if line.find("Molecular mass for molecule type     1 is") != -1:
                        molecMass1 = line.replace("Molecular mass for molecule type     1 is", "")
                        molecMass1 = molecMass1.replace("g/mol", "")
                        molecMass1 = molecMass1.strip()
                        molecMassList1.append(float(molecMass1))

                    if line.find("Molecule Number                1") != -1:
                        molecNum = line.replace("Molecule Number                1", "")
                        molecNum = molecNum.strip()
                        molecNumBox1.append(float(molecNum[0:8]))
                        molecNumBox2.append(float(molecNum[11:19]))

                    if line.find("Molar Volume         ml/mol") != -1:
                        molarVolume = line.replace("Molar Volume         ml/mol", "")
                        molarVolume = molarVolume.strip()
                        molarVolumeBox1.append(float(molarVolume[0:11]))
                        molarVolumeBox2.append(float(molarVolume[12:23]))

                    # Pressure is converted from kPa back to bar.

                    if line.find("Virial Pressure         kPa") != -1:
                        virialPress = line.replace("Virial Pressure         kPa", "")
                        virialPress = virialPress.strip()
                        virialPressBox1.append(float(virialPress[0:11]) / 100)
                        virialPressBox2.append(float(virialPress[12:23]) / 100)

                    # Pressure is converted from kPa back to bar.

                    if line.find("Thermodynamic Pressure  kPa") != -1:
                        thermoPress = line.replace("Thermodynamic Pressure  kPa", "")
                        thermoPress = thermoPress.strip()
                        thermoPressBox1.append(float(thermoPress[0:11]) / 100)
                        thermoPressBox2.append(float(thermoPress[12:23]) / 100)

                    if line.find("u (Gibbs Total)           K    1") != -1:
                        chemPot = line.replace("u (Gibbs Total)           K    1", "")
                        chemPot = chemPot.strip()
                        chemPotBox1.append(float(chemPot[0:11]))
                        chemPotBox2.append(float(chemPot[12:23]))

                    if line.find("H_vap (vapor p)      kJ/mol") != -1:
                        deltaH = line.replace("H_vap (vapor p)      kJ/mol", "")
                        deltaH = deltaH.strip()
                        deltaHList.append(float(deltaH))

        # The 'currentPath' from earlier is used to open an output file in the initial directory for printing select
        # numerical results and averages. Note that statistical computations are performed while the file is open,
        # and while this is not ideal for efficiency or habit it's trivial for this code.

        os.chdir(currentPath)
        fileOutput = open("towhee_" + speciesDir, "w")

        # Depending on the name of the output file, equilibration or production is written to the output file.

        fileOutput.write("\n")
        if towheeOutput == "towhee_output1.txt":
            fileOutput.write("** Equilibration Results **\n")

        elif towheeOutput == "towhee_output2.txt":
            fileOutput.write("** Production Results **\n")

        # The name of the output file (each file has the same file name) is written to the output file.

        fileOutput.write("\n")
        fileOutput.write("File: " + towheeOutput + "\n")
        fileOutput.write("\n")

        # The ensemble type is written to the file.

        fileOutput.write("Ensemble: " + ensemble + "\n")
        fileOutput.write("\n")

        # The system name is written to the output file via 'speciesDir'.

        fileOutput.write("System Type: " + speciesDir + "\n")
        fileOutput.write("\n")

        # The average number of types of molecules in each set it written to the file.

        avgNumMolTy = int(statistics.mean(nMolTyList))
        fileOutput.write("Number of Molecule Type: " + str(avgNumMolTy) + "\n")
        fileOutput.write("\n")

        # A string specifying location of standard deviations in written to the file.

        fileOutput.write("** Standard deviations are shown in parentheses **\n")
        fileOutput.write("\n")

        # The number of sets is written to the file.

        fileOutput.write("Number of Sets:\n")
        fileOutput.write(str(len(setList)))
        fileOutput.write("\n")

        # The 'setList' is written to the file.

        fileOutput.write("\n")
        fileOutput.write("Sets:\n")
        fileOutput.write(str(setList))
        fileOutput.write("\n")

        # The 'tempList' is written to the file.

        fileOutput.write("\n")
        fileOutput.write("Set Temperature in K:\n")
        fileOutput.write(str(temperatureList))
        fileOutput.write("\n")

        if ensemble != "gibbs":

            # The 'pressureList' is written to the file.

            fileOutput.write("\n")
            fileOutput.write("Set Pressure in bar:\n")
            fileOutput.write(str(pressureList))
            fileOutput.write("\n")

            # The 'nMolecList' is written to the file.

            fileOutput.write("\n")
            fileOutput.write("Set Total Number of Molecules:\n")
            fileOutput.write(str(nMolecList))
            fileOutput.write("\n")

            # The average number of molecules is written to the file.

            fileOutput.write("\n")
            avgNumMolec = statistics.mean(nMolecList)
            fileOutput.write("Average Total Number of Molecules:\n")
            fileOutput.write(str(int(avgNumMolec)))
            fileOutput.write("\n")

            # 'nMolecTypList' is written to the file. This is a self consistency check and should always be equal to
            # the average number of molecules for pure components.

            fileOutput.write("\n")
            fileOutput.write("Number of Molecules of Each Type:\n")
            fileOutput.write(str(nMolecTypList))
            fileOutput.write("\n")

            # The molar mass of component 1 is written to the file.

            fileOutput.write("\n")
            fileOutput.write("Molecular Mass for Molecule Type 1 in g/mol: \n")
            fileOutput.write(str(molecMassList1))
            fileOutput.write("\n")

            # If the system is a binary, the molar mass of component 2 is written to the file.

            if avgNumMolTy > 1:
                fileOutput.write("\n")
                fileOutput.write("Molecular Mass for Molecule Type 2 in g/mol: \n")
                fileOutput.write(str(molecMassList2))
                fileOutput.write("\n")

            # Internal energy of departure list 'internalEnergyList' and the average in various units is
            # written to the file.

            fileOutput.write("\n")
            fileOutput.write("System Internal Energies of Departure in kJ/mol:\n")
            fileOutput.write(str(internalEnergyList))
            fileOutput.write("\n")

            avgInternalEnergy = statistics.mean(internalEnergyList)
            stdevInternalEnergy = statistics.stdev(internalEnergyList)
            pStdevInternalEnergy = abs(stdevInternalEnergy / avgInternalEnergy * 100.0)
            fileOutput.write("\n")
            fileOutput.write("System Internal Energy of Departure in kJ/mol:\n")
            fileOutput.write(str(avgInternalEnergy) + "(+/- " + str(stdevInternalEnergy) + ")" +
                             "(+/- " + str(pStdevInternalEnergy) + "%)\n")
            fileOutput.write("\n")

            avgInternalEnergy *= 10000
            stdevInternalEnergy *= 10000
            pStdevInternalEnergy = abs(stdevInternalEnergy / avgInternalEnergy * 100.0)
            fileOutput.write("System Internal Energy of Departure in ccbar/mol:\n")
            fileOutput.write(str(avgInternalEnergy) + "(+/- " + str(stdevInternalEnergy) + ")" +
                             "(+/- " + str(pStdevInternalEnergy) + "%)\n")
            fileOutput.write("\n")

            avgInternalEnergy = avgInternalEnergy / avgNumMolec
            stdevInternalEnergy = stdevInternalEnergy / avgNumMolec
            pStdevInternalEnergy = abs(stdevInternalEnergy / avgInternalEnergy * 100.0)
            fileOutput.write("Internal Energy of Departure (U/N) in ccbar/mol:\n")
            fileOutput.write(str(avgInternalEnergy) + "(+/- " + str(stdevInternalEnergy) + ")" +
                             "(+/- " + str(pStdevInternalEnergy) + "%)\n")
            fileOutput.write("\n")

            # The list of set molar volumes 'molarVolumeList' and the average molar volume is
            # written to the file

            fileOutput.write("Average Molar Volumes in cc/mol:\n")
            fileOutput.write(str(molarVolumeList))
            fileOutput.write("\n")

            avgMolarVolume = statistics.mean(molarVolumeList)
            stdevMolarVolume = statistics.stdev(molarVolumeList)
            pStdevMolarVolume = abs(stdevMolarVolume / avgMolarVolume * 100.0)
            fileOutput.write("\n")
            fileOutput.write("Average Molar Volume in cc/mol:\n")
            fileOutput.write(str(avgMolarVolume) + "(+/- " + str(stdevMolarVolume) + ")" +
                             "(+/- " + str(pStdevMolarVolume) + "%)\n")

            # For single component systems, the average density in various units is written to the file.

            if avgNumMolTy == 1:
                avgMolecMass1 = statistics.mean(molecMassList1)
                avgMolarDensity = 1.0 / avgMolarVolume
                avgMassDensity1 = avgMolarDensity * avgMolecMass1
                avgMassDensity2 = avgMassDensity1 * 1000.0

                fileOutput.write("\n")
                fileOutput.write("Average Molar Density in mol/cc: \n")
                fileOutput.write(str(avgMolarDensity))
                fileOutput.write("\n")

                fileOutput.write("\n")
                fileOutput.write("Average Mass Density in g/cc: \n")
                fileOutput.write(str(avgMassDensity1))
                fileOutput.write("\n")

                fileOutput.write("\n")
                fileOutput.write("Average Mass Density in kg/m3: \n")
                fileOutput.write(str(avgMassDensity2))
                fileOutput.write("\n")

        # A flag is written to the file if the average number of molecules is not equal to the number
        # of molecules in each set. This is the only flag currently in the code.

            for eachMolecNum in nMolecList:
                if eachMolecNum != avgNumMolec:
                    fileOutput.write("** Number of Molecules are not equal among sets **\n")

        if ensemble == "gibbs":

            # The 'nMolecList' is written to the file.

            fileOutput.write("\n")
            fileOutput.write("Set Total Number of Molecules:\n")
            fileOutput.write(str(nMolecList))
            fileOutput.write("\n")

            # The average number of molecules is written to the file.

            fileOutput.write("\n")
            avgNumMolec = statistics.mean(nMolecList)
            fileOutput.write("Average Total Number of Molecules:\n")
            fileOutput.write(str(avgNumMolec))
            fileOutput.write("\n")

            # The molar mass of component 1 is written to the file.

            fileOutput.write("\n")
            fileOutput.write("Molecular Mass for Molecule Type 1 in g/mol: \n")
            fileOutput.write(str(molecMassList1))
            fileOutput.write("\n")

            # The number of molecules in box 1 is written to the file.

            fileOutput.write("\n")
            fileOutput.write("Number of Molecules Box 1: \n")
            fileOutput.write(str(molecNumBox1))
            fileOutput.write("\n")

            # The average number of molecules in box 1 is written to the file

            fileOutput.write("\n")
            avgNumMolecBox1 = statistics.mean(molecNumBox1)
            stdevNumMolecBox1 = statistics.stdev(molecNumBox1)
            pStdevNumMolecBox1 = abs(stdevNumMolecBox1 / avgNumMolecBox1 * 100.0)
            fileOutput.write("Average Number of Molecules Box 1:\n")
            fileOutput.write(str(avgNumMolecBox1) + "(+/- " + str(stdevNumMolecBox1) + ")" +
                             "(+/- " + str(pStdevNumMolecBox1) + "%)")
            fileOutput.write("\n")

            # The number of molecules in box 2 is written to the file.

            fileOutput.write("\n")
            fileOutput.write("Number of Molecules Box 2: \n")
            fileOutput.write(str(molecNumBox2))
            fileOutput.write("\n")

            # The average number of molecules in box 2 is written to the file

            fileOutput.write("\n")
            avgNumMolecBox2 = statistics.mean(molecNumBox2)
            stdevNumMolecBox2 = statistics.stdev(molecNumBox2)
            pStdevNumMolecBox2 = abs(stdevNumMolecBox2 / avgNumMolecBox2 * 100.0)
            fileOutput.write("Average Number of Molecules Box 2:\n")
            fileOutput.write(str(avgNumMolecBox2) + "(+/- " + str(stdevNumMolecBox2) + ")" +
                             "(+/- " + str(pStdevNumMolecBox2) + "%)")
            fileOutput.write("\n")

            # The molar volume in box 1 is written to the file.

            fileOutput.write("\n")
            fileOutput.write("Molar Volume Box 1 in cc/mol: \n")
            fileOutput.write(str(molarVolumeBox1))
            fileOutput.write("\n")

            # The average molar volume in box 1 is written to the file

            fileOutput.write("\n")
            avgMolarVolumeBox1 = statistics.mean(molarVolumeBox1)
            stdevMolarVolumeBox1 = statistics.stdev(molarVolumeBox1)
            pStdevMolarVolumeBox1 = abs(stdevMolarVolumeBox1 / avgMolarVolumeBox1 * 100.0)
            fileOutput.write("Average Molar Volume Box 1 in cc/mol:\n")
            fileOutput.write(str(avgMolarVolumeBox1) + "(+/- " + str(stdevMolarVolumeBox1) + ")" +
                             "(+/- " + str(pStdevMolarVolumeBox1) + "%)")
            fileOutput.write("\n")

            # The molar volume in box 2 is written to the file.

            fileOutput.write("\n")
            fileOutput.write("Molar Volume Box 2 in cc/mol: \n")
            fileOutput.write(str(molarVolumeBox2))
            fileOutput.write("\n")

            # The average molar volume in box 2 is written to the file

            fileOutput.write("\n")
            avgMolarVolumeBox2 = statistics.mean(molarVolumeBox2)
            stdevMolarVolumeBox2 = statistics.stdev(molarVolumeBox2)
            pStdevMolarVolumeBox2 = abs(stdevMolarVolumeBox2 / avgMolarVolumeBox2 * 100.0)
            fileOutput.write("Average Molar Volume Box 2 in cc/mol:\n")
            fileOutput.write(str(avgMolarVolumeBox2) + "(+/- " + str(stdevMolarVolumeBox2) + ")" +
                             "(+/- " + str(pStdevMolarVolumeBox2) + "%)")
            fileOutput.write("\n")

            # The virial pressure in box 1 is written to the file.

            fileOutput.write("\n")
            fileOutput.write("Virial Pressure Box 1 in Bar: \n")
            fileOutput.write(str(virialPressBox1))
            fileOutput.write("\n")

            # The average virial pressure in box 1 is written to the file

            fileOutput.write("\n")
            avgVirialPressBox1 = statistics.mean(virialPressBox1)
            stdevVirialPressBox1 = statistics.stdev(virialPressBox1)
            pStdevVirialPressBox1 = abs(stdevVirialPressBox1 / avgVirialPressBox1 * 100.0)
            fileOutput.write("Average Virial Pressure Box 1 in Bar:\n")
            fileOutput.write(str(avgVirialPressBox1) + "(+/- " + str(stdevVirialPressBox1) + ")" +
                             "(+/- " + str(pStdevVirialPressBox1) + "%)")
            fileOutput.write("\n")

            # The virial pressure in box 2 is written to the file.

            fileOutput.write("\n")
            fileOutput.write("Virial Pressure Box 2 in Bar: \n")
            fileOutput.write(str(virialPressBox2))
            fileOutput.write("\n")

            # The average virial pressure in box 2 is written to the file

            fileOutput.write("\n")
            avgVirialPressBox2 = statistics.mean(virialPressBox2)
            stdevVirialPressBox2 = statistics.stdev(virialPressBox2)
            pStdevVirialPressBox2 = abs(stdevVirialPressBox2 / avgVirialPressBox2 * 100.0)
            fileOutput.write("Average Virial Pressure Box 2 in Bar:\n")
            fileOutput.write(str(avgVirialPressBox2) + "(+/- " + str(stdevVirialPressBox2) + ")" +
                             "(+/- " + str(pStdevVirialPressBox2) + "%)")
            fileOutput.write("\n")

            # The thermodynamic pressure in box 1 is written to the file.

            fileOutput.write("\n")
            fileOutput.write("Thermodynamic Pressure Box 1 in Bar: \n")
            fileOutput.write(str(thermoPressBox1))
            fileOutput.write("\n")

            # The average thermodynamic pressure in box 1 is written to the file

            fileOutput.write("\n")
            avgThermoPressBox1 = statistics.mean(thermoPressBox1)
            stdevThermoPressBox1 = statistics.stdev(thermoPressBox1)
            pStdevThermoPressBox1 = abs(stdevThermoPressBox1 / avgThermoPressBox1 * 100.0)
            fileOutput.write("Average Thermodynamic Pressure Box 1 in Bar:\n")
            fileOutput.write(str(avgThermoPressBox1) + "(+/- " + str(stdevThermoPressBox1) + ")" +
                             "(+/- " + str(pStdevThermoPressBox1) + "%)")
            fileOutput.write("\n")

            # The thermodynamic pressure in box 2 is written to the file.

            fileOutput.write("\n")
            fileOutput.write("Thermodynamic Pressure Box 2 in Bar: \n")
            fileOutput.write(str(thermoPressBox2))
            fileOutput.write("\n")

            # The average thermodynamic pressure in box 2 is written to the file

            fileOutput.write("\n")
            avgThermoPressBox2 = statistics.mean(thermoPressBox2)
            stdevThermoPressBox2 = statistics.stdev(thermoPressBox2)
            pStdevThermoPressBox2 = abs(stdevThermoPressBox2 / avgThermoPressBox2 * 100.0)
            fileOutput.write("Average Thermodynamic Pressure Box 2 in Bar:\n")
            fileOutput.write(str(avgThermoPressBox2) + "(+/- " + str(stdevThermoPressBox2) + ")" +
                             "(+/- " + str(pStdevThermoPressBox2) + "%)")
            fileOutput.write("\n")

            # The total gibbs chemical potential in box 1 is written to the file.

            fileOutput.write("\n")
            fileOutput.write("Total Gibbs Chemical Potential Box 1 in K: \n")
            fileOutput.write(str(chemPotBox1))
            fileOutput.write("\n")

            # The average total gibbs chemical potential in box 1 is written to the file

            fileOutput.write("\n")
            avgChemPotBox1 = statistics.mean(chemPotBox1)
            stdevChemPotBox1 = statistics.stdev(chemPotBox1)
            pStdevChemPotBox1 = abs(stdevChemPotBox1 / avgChemPotBox1 * 100.0)
            fileOutput.write("Average Total Gibbs Chemical Potential Box 1 in K:\n")
            fileOutput.write(str(avgChemPotBox1) + "(+/- " + str(stdevChemPotBox1) + ")" +
                             "(+/- " + str(pStdevChemPotBox1) + "%)")
            fileOutput.write("\n")

            # The chemical potential in box 2 is written to the file.

            fileOutput.write("\n")
            fileOutput.write("Total Gibbs Chemical Potential Box 2 in K: \n")
            fileOutput.write(str(chemPotBox2))
            fileOutput.write("\n")

            # The average total gibbs chemical potential in box 1 is written to the file

            fileOutput.write("\n")
            avgChemPotBox2 = statistics.mean(chemPotBox2)
            stdevChemPotBox2 = statistics.stdev(chemPotBox2)
            pStdevChemPotBox2 = abs(stdevChemPotBox2 / avgChemPotBox2 * 100.0)
            fileOutput.write("Average Total Gibbs Chemical Potential Box 2 in K:\n")
            fileOutput.write(str(avgChemPotBox2) + "(+/- " + str(stdevChemPotBox2) + ")" +
                             "(+/- " + str(pStdevChemPotBox2) + "%)")
            fileOutput.write("\n")

            # The chemical potential percent difference is written to the file

            fileOutput.write("\n")
            percentDifChemPot = 2 * abs(avgChemPotBox1 - avgChemPotBox2) / abs(avgChemPotBox1 + avgChemPotBox2) * 100
            fileOutput.write("Chemical Potential Percent Difference: \n")
            fileOutput.write(str(percentDifChemPot) + "%")
            fileOutput.write("\n")

            # The enthalpy of vaporization is written to the file.

            fileOutput.write("\n")
            fileOutput.write("System Enthalpy of Vaporization in kJ/mol: \n")
            fileOutput.write(str(deltaHList))
            fileOutput.write("\n")

            # The average enthalpy of vaporization is written to the file

            fileOutput.write("\n")
            avgDeltaH = statistics.mean(deltaHList)
            stdevDeltaH = statistics.stdev(deltaHList)
            pStdevDeltaH = abs(stdevDeltaH / avgDeltaH * 100.0)
            fileOutput.write("Average System Enthalpy of Vaporization in kJ/mol:\n")
            fileOutput.write(str(avgDeltaH) + "(+/- " + str(stdevDeltaH) + ")" +
                             "(+/- " + str(pStdevDeltaH) + "%)")
            fileOutput.write("\n")

            # A flag is written to the file if the average number of molecules is not equal to the number
            # of molecules in each set. This is the only flag currently in the code.

            for eachMolecNum in nMolecList:
                if eachMolecNum != avgNumMolec:
                    fileOutput.write("** Number of Molecules are not equal among sets **\n")

        # The output file is closed.

        fileOutput.close()

# End towhee_collect
