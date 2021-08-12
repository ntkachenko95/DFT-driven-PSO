The following python libraries are required: numpy, scipy, matplotlib
The code is compatible with the Gaussian16 code. Compatibility with other software will be added in the future.
Description of important files and folders:

1) /scripts/Launcher.py - main script that needs to be launched by a cluster computing node
2) /scripts/Script_creation.py - File that creates Slurm scripts and should be adjusted for your needs. Each HPC-center uses its own sbatch scripts. Please replace mine with the appropriate one.
3) /Input_files folder contains all necessary input files for the test run calculation. Extra solvent molecules might be specified if needed.
4) /Input_files/Input.txt - main input file with important parameters of DFT-driven-PSO. The provided Input.txt file is as an example, and the parameters might not be optimal.

To launch the program, copy all scripts and input files into one folder. Do not forget to change the "main_path" in the Input.txt file. 
Run Launcher.py on the computing node with the command "python Launcher.py".
