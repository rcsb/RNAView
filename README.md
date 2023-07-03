# RNAView
RNAView identifies base pairs that are formed in nucleic acid structures and classifies them according to the system of Leontis and Westhof

## Introduction
RNAVIEW program quickly displays the secondary structure of RNA/DNA with tertiary interactions. It fully implements Leontis and Westhof's convention for the edge-to-edge hydrogen bonding interactions. 

Leontis NB and Westhof E. (2001) Geometric nomenclature and classification of RNA base pairs. RNA 7:499-512. https://doi.org/10.1017/s1355838201002515

The software has been tested in Linux and UNIX (SGI), MAC, SUN systems

## Installation
You can install the program in any directory.

1. unpack the program package in the desired directory (e.g. /?/?/?/ ).

	```
	zcat  RNAVIEW.tar.gz | tar xvf -
 	```

2. install the program
   
   ```
   cd /?/?/?/RNAVIEW/
   make
   ```
   
   The executable file rnaview will be in the directory  /?/?/?/RNAVIEW/bin/.

3. Define RNAVIEW environment variable to point to the installation directory (e.g. /?/?/?/RNAVIEW)
   Add the following command to your shell script

   For C shell users:

   ```
    setenv RNAVIEW /?/?/?/RNAVIEW
    setenv PATH "/?/?/?/RNAVIEW/bin:"$PATH
    ```

   For Bourne shell users:

   ```
   RNAVIEW=/?/?/?/RNA/RNAVIEW; export RNAVIEW
   PATH="/?/?/?/RNAVIEW/bin:"$PATH; export PATH
   ```

4. To test the program, go to /?/?/?/RNAVIEW/test

   ```
   rnaview -p tr0001.pdb
   ```
   
   You should get a postscript file *.ps and some other output files.



## Usage

type rnaview  or  rnaview -h  to  get an online help.

Usage: executable [option]  input_file

## Options

1. If no [option] is given, it will only generate the fully annotated base pair lists. Example:    rnaview  pdbfile_name
2. [option] -p to generate fully annotated 2D structure in postscript format. Detailed information is given in XML format(RNAML) Example: rnaview  -p pdbfile_name                                                 
3. [option] -v to generate a 3D structure in VRML format.  It can be displayed on internet (with VRML plug in). Example:    rnaview  -v pdbfile_name
4. [option] -c to select chains for calculation. -c should be followed by chain IDs. If select several chains, they should be put together, like ABC for chain A, B and C.  This option is useful, when drawing a single copy of 2D structure from a dimer or trimer PDB file.  Example:    rnaview  -pc ABC pdbfile_name      
5. [option] -a to process many pdbfiles. The pdbfile names must be put in one file (like  file.list) and seperated by a space. You may give the resolution after file.list. If you do not give (or give 0), it means resolution is  ignored   Example:    rnaview  -a file.list 3.0   means that only the pdbfiles with resolution < 3.0  are selected for calculation.
6. [option] -x to input XML (RNAML) file.  Normally this option is combined with -p to generate a 2D structure.  Example:    rnaview  -px RNAML_file_name

## NMR structures

The program will automatically pick the best NMR model according to the REMARK record of the pdb file. If there is no best model in the PDB file, it will pick the first model from the ensemble.


