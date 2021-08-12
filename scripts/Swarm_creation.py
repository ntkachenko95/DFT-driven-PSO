import numpy as np
import operator
import shutil
import os
from scipy.optimize import fsolve
import copy
import collections
import matplotlib.pyplot as plt
from Script_creation import *
from PSO_tools import *



#######################Define all variables#####################################

#reading variables from the Input.txt file

Generation=1

with open('Input.txt', 'r') as f:
    file = f.read()
file = file.split('\n')
file = [i.split() for i in file]

main_path = file[2][1]
script_file_header = file[3][1]
num_generations = int(file[4][1])

Generation_size = int(file[8][1])
Saved_population = int(file[9][1])
Box_size = float(file[10][1])
Is_Solvatation = bool(file[11][1])
fragment_separation = float(file[12][1])
fragment_types = int(file[13][1])
fragment_numbers = list(map(int,file[14][1:]))
fragment_filenames = []
for i in range(fragment_types):
    fragment_filenames.append(file[15+i][1])
Header_for_gaussian_DFT_filename = file[15+fragment_types][1]
Header_for_gaussian_PM3_filename = file[16+fragment_types][1]


#additional readings

with open(Header_for_gaussian_DFT_filename, 'r') as f:
    Header_for_gaussian = f.read()
with open(Header_for_gaussian_PM3_filename, 'r') as f:
    Header_for_gaussian_PM3 = f.read()
fragments = []
for filenames_xyz in fragment_filenames:
    with open(filenames_xyz, 'r') as f:
        xyz_string = f.read()
    fragments.append(fragment_alignment(fragment_from_string(xyz_string)))


#Create_swarm
Swarm=[]
for i in range(Generation_size):
    Swarm.append(Particle('particle_'+str(i), fragment_numbers, fragments, Box_size=Box_size, Is_Solvatation=Is_Solvatation))
    Swarm[i].random_initialization()
    Swarm[i].fix_overlap()
    make_gaussian_file(Swarm[i], Header_for_gaussian_PM3)

for i in range(Generation_size):
    Script_creation(script_file_header+'_particle_'+str(i), i, main_path, Generation)

os.mkdir('Generation_'+str(Generation)+'_P3')
for i in range(Generation_size):
    shutil.copyfile(script_file_header+'_particle_'+str(i), 'Generation_'+str(Generation)+'_P3/'+script_file_header+'_particle_'+str(i))
    shutil.copyfile('particle_'+str(i)+'.com', 'Generation_'+str(Generation)+'_P3/particle_'+str(i)+'.com')

import pickle
with open('Swarm.data', 'wb') as f:
    pickle.dump(Swarm, f)

with open('Generation.data', 'wb') as f:
    pickle.dump(Generation, f)