import numpy as np
import operator
import shutil
import os
from scipy.optimize import fsolve
import copy
import collections
import matplotlib.pyplot as plt
from Script_creation import *

def make_gaussian_file_from_preoptimization(filename,Header, Coords):
    with open(filename, 'w') as f:
        f.write(Header+'\n')
        f.write(Coords+'\n')
def fragment_from_string(string):
    fragment = string.split('\n')
    for i in range(len(fragment)):
        fragment[i]=fragment[i].split()
        fragment[i][0]=int(fragment[i][0])
        fragment[i][1]=float(fragment[i][1])
        fragment[i][2]=float(fragment[i][2])
        fragment[i][3]=float(fragment[i][3])
    return fragment
def string_from_fragment(fragment):
    a=str()
    for i in fragment:
        a+=str(int(i[0]))+' '+str(i[1])+' '+str(i[2])+' '+str(i[3])+'\n'
    return(a)
#fragment - list [[Atomic_number, x, y, z], ...]
def fragment_center_of_mass(fragment):
    X0, Y0, Z0 = (sum([n[0]*n[1] for n in fragment])/sum([n[0] for n in fragment]),sum([n[0]*n[2] for n in fragment])/sum([n[0] for n in fragment]),sum([n[0]*n[3] for n in fragment])/sum([n[0] for n in fragment]))
    return [X0, Y0, Z0]

def fragment_alignment(fragment):
    X0, Y0, Z0 = (sum([n[0]*n[1] for n in fragment])/sum([n[0] for n in fragment]),sum([n[0]*n[2] for n in fragment])/sum([n[0] for n in fragment]),sum([n[0]*n[3] for n in fragment])/sum([n[0] for n in fragment]))
    new_fragment = copy.deepcopy(fragment)
    for i in range(len(new_fragment)):
        new_fragment[i][1]-=X0
        new_fragment[i][2]-=Y0
        new_fragment[i][3]-=Z0
    return new_fragment        
def fragment_rotation_xyz(fragment,alpha,beta,gamma):
    new_fragment = copy.deepcopy(fragment)
    for i in range(len(new_fragment)):
        (new_fragment[i][1],new_fragment[i][2],new_fragment[i][3])=(np.cos(alpha)*np.cos(beta)*new_fragment[i][1]+(np.cos(alpha)*np.sin(beta)*np.sin(gamma)-np.sin(alpha)*np.cos(gamma))*new_fragment[i][2]+(np.cos(alpha)*np.sin(beta)*np.cos(gamma)+np.sin(alpha)*np.sin(gamma))*new_fragment[i][3],np.sin(alpha)*np.cos(beta)*new_fragment[i][1]+(np.sin(alpha)*np.sin(beta)*np.sin(gamma)+np.cos(alpha)*np.cos(gamma))*new_fragment[i][2]+(np.sin(alpha)*np.sin(beta)*np.cos(gamma)-np.cos(alpha)*np.sin(gamma))*new_fragment[i][3],-np.sin(beta)*new_fragment[i][1]+np.cos(beta)*np.sin(gamma)*new_fragment[i][2]+np.cos(beta)*np.cos(gamma)*new_fragment[i][3])
    return new_fragment

def fragment_shift(fragment, x, y, z):
    new_fragment = copy.deepcopy(fragment)
    for i in range(len(new_fragment)):
        new_fragment[i][1]-=x
        new_fragment[i][2]-=y
        new_fragment[i][3]-=z
    return new_fragment
def fragment_placement(fragment, x ,y ,z):
    new_fragment = copy.deepcopy(fragment)
    for i in range(len(new_fragment)):
        new_fragment[i][1]+=x
        new_fragment[i][2]+=y
        new_fragment[i][3]+=z
    return new_fragment
def fragment_radii(fragment):
    new_fragment=fragment_alignment(copy.deepcopy(fragment))
    radii=0
    for i in copy.deepcopy(fragment):
        dist=np.sqrt(i[1]**2+i[2]**2+i[3]**2)
        if dist>=radii: radii=dist
    return(dist)
def freeze_molecule(number_of_atoms, number_of_atoms_accounted):
    freeze_string=str()
    for i in range(number_of_atoms-1):
        for j in range(i+1,number_of_atoms):
            freeze_string+='R('+str(i+1+number_of_atoms_accounted)+','+str(j+1+number_of_atoms_accounted)+') freeze \n'
    number_of_atoms_new=number_of_atoms_accounted+number_of_atoms
    return freeze_string, number_of_atoms_new

def slow_molecule(number_of_atoms, number_of_atoms_accounted):
    freeze_string=str()
    for i in range(1,number_of_atoms):
        freeze_string+='R('+str(1+number_of_atoms_accounted)+','+str(i+1+number_of_atoms_accounted)+') freeze \n'
    number_of_atoms_new=number_of_atoms_accounted+number_of_atoms
    return freeze_string, number_of_atoms_new

def make_gaussian_file(particle, Header):
    with open(particle.name+'.com', 'w') as f:
        f.write(Header+'\n')
        f.write(particle.print_coordinates()+'\n')
        #number_of_atoms_accounted=0
        #for specie_index in range(particle.num_species):
        #    for fragment_number in range(particle.fragment_numbers[specie_index]):
        #        freeze_string, number_of_atoms_accounted = slow_molecule(len(particle.fragments[specie_index]), number_of_atoms_accounted)
        #        f.write(freeze_string)
def RMS_distance_evaluation(particle_1,particle_2):
    coord1=copy.deepcopy(particle_1.get_full_coordinates())
    coord2=copy.deepcopy(particle_2.get_full_coordinates())
    number_of_atoms=len(coord1)
    CM_1=fragment_center_of_mass(coord1)
    CM_2=fragment_center_of_mass(coord2)
    for i in range(len(coord1)):
        coord1[i].append(np.sqrt((coord1[i][1])**2+(coord1[i][2])**2+(coord1[i][3])**2))
    coord1 = sorted(coord1, key=operator.itemgetter(0, 4))
    for i in range(len(coord2)):
        coord2[i].append(np.sqrt((coord2[i][1])**2+(coord2[i][2])**2+(coord2[i][3])**2))
    coord2 = sorted(coord2, key=operator.itemgetter(0, 4))
    D_I=np.zeros((number_of_atoms,number_of_atoms))
    D_II=np.zeros((number_of_atoms,number_of_atoms))
    for i in range(number_of_atoms-1):
        for j in range(i+1,number_of_atoms):
            D_I[i][j]=np.sqrt((coord1[i][1]-coord1[j][1])**2+(coord1[i][2]-coord1[j][2])**2+(coord1[i][3]-coord1[j][3])**2)
            D_II[i][j]=np.sqrt((coord2[i][1]-coord2[j][1])**2+(coord2[i][2]-coord2[j][2])**2+(coord2[i][3]-coord2[j][3])**2)
    Square=0
    for i in range(number_of_atoms-1):
        for j in range(i,number_of_atoms):
            Square+=(D_I[i][j]-D_II[i][j])**2
    RMS = np.sqrt((2/(number_of_atoms**2-number_of_atoms))*Square)
    return RMS

def other_RMS_distance_evaluation(particle_1,particle_2):
    coord1=copy.deepcopy(particle_1.get_full_coordinates())
    coord2=copy.deepcopy(particle_2.get_full_coordinates())
    number_of_atoms=len(coord1)
    CM_1=fragment_center_of_mass(coord1)
    CM_2=fragment_center_of_mass(coord2)
    for i in range(len(coord1)):
        coord1[i].append(np.sqrt((coord1[i][1])**2+(coord1[i][2])**2+(coord1[i][3])**2))
    coord1 = sorted(coord1, key=operator.itemgetter(0, 4))
    for i in range(len(coord2)):
        coord2[i].append(np.sqrt((coord2[i][1])**2+(coord2[i][2])**2+(coord2[i][3])**2))
    coord2 = sorted(coord2, key=operator.itemgetter(0, 4))
    Square=0
    for i in range(number_of_atoms):
        Square+=(coord1[i][4]-coord2[i][4])**2
    RMS = np.sqrt((1/(number_of_atoms))*Square)
    return RMS
def angles_from_coordinates(x0,y0,z0,x1,y1,z1):
    def func(x):
        return [np.cos(x[0])*np.cos(x[1])*x0+(np.cos(x[0])*np.sin(x[1])*np.sin(x[2])-np.sin(x[0])*np.cos(x[2]))*y0+(np.cos(x[0])*np.sin(x[1])*np.cos(x[2])+np.sin(x[0])*np.sin(x[2]))*z0-x1,
               np.sin(x[0])*np.cos(x[1])*x0+(np.sin(x[0])*np.sin(x[1])*np.sin(x[2])+np.cos(x[0])*np.cos(x[2]))*y0+(np.sin(x[0])*np.sin(x[1])*np.cos(x[2])-np.cos(x[0])*np.sin(x[2]))*z0-y1,
               -np.sin(x[1])*x0+np.cos(x[1])*np.sin(x[2])*y0+np.cos(x[1])*np.cos(x[2])*z0 - z1]

    root = fsolve(func, [0,0,0])
    alpha, beta, gamma=root
    return(np.array([alpha, beta, gamma]))

def average_angle_from_fragment(initial_fragment, fragment):
    angles = []
    for index,atom in enumerate(copy.deepcopy(fragment)):
        x0=atom[1]
        y0=atom[2]
        z0=atom[3]
        if np.sqrt(x0**2+y0**2+z0**2)<0.1:
            pass
        else:
            angles.append(angles_from_coordinates(initial_fragment[index][1],initial_fragment[index][2],initial_fragment[index][3],atom[1],atom[2],atom[3]))
    for i in angles:
        i=[i[0]%(2*np.pi),i[1]%(2*np.pi),i[2]%(2*np.pi)]
    
    return [sum([n[0] for n in angles])/len(angles),sum([n[1] for n in angles])/len(angles),sum([n[2] for n in angles])/len(angles)]

def angles_from_fragment(initial_fragment,fragment):
    def Large_func(x):
        result=[]
        for index,atom in enumerate(copy.deepcopy(fragment)):
            x0=initial_fragment[index][1]
            y0=initial_fragment[index][2]
            z0=initial_fragment[index][3]
            x1=atom[1]
            y1=atom[2]
            z1=atom[3]
            if np.sqrt(x0**2+y0**2+z0**2)<0.1:
                pass
            else:
                result.append([(np.cos(x[0])*np.cos(x[1])*x0+(np.cos(x[0])*np.sin(x[1])*np.sin(x[2])-np.sin(x[0])*np.cos(x[2]))*y0+(np.cos(x[0])*np.sin(x[1])*np.cos(x[2])+np.sin(x[0])*np.sin(x[2]))*z0-x1)**2,
                   (np.sin(x[0])*np.cos(x[1])*x0+(np.sin(x[0])*np.sin(x[1])*np.sin(x[2])+np.cos(x[0])*np.cos(x[2]))*y0+(np.sin(x[0])*np.sin(x[1])*np.cos(x[2])-np.cos(x[0])*np.sin(x[2]))*z0-y1)**2,
                   (-np.sin(x[1])*x0+np.cos(x[1])*np.sin(x[2])*y0+np.cos(x[1])*np.cos(x[2])*z0 - z1)**2])
        return([sum([n[0] for n in result]),sum([n[1] for n in result]),sum([n[2] for n in result])])
    root = fsolve(Large_func, [0,0,0], factor=0.5)
    alpha, beta, gamma=[root[0]%(2*np.pi),root[1]%(2*np.pi),root[2]%(2*np.pi)]
    #print(np.isclose(Large_func(root), [0.0, 0.0, 0.0]))
    #print(root)
    #print(Large_func(root))
    #print(Large_func([root[0]%(2*np.pi),root[1]%(2*np.pi),root[2]%(2*np.pi)]))
    #print('***************')
    return(np.array([alpha, beta, gamma]))

def read_gaussian_file(particle):
    Normal_termination = True
    optimized_coordinates = []
    optimized_particle_coordinates = []
    with open(particle.name+'.log', 'r') as f:
        file = f.read()
        file = file.split('\n')
        if 'Error' in file[-5]:
            Normal_termination = False
        coordinates_index=len(file)-1-file[::-1].index('                          Input orientation:                          ')
        energy_index = -1
        for i in file[::-1]:
            energy_index+=1
            if 'SCF Done:' in i:
                energy_index=len(file)-1-energy_index
                break
        energy = float(file[energy_index].split()[4])
        counter=0
        for specie_index in range(particle.num_species):
            optimized_coordinates.append([])
            for fragment_number in range(particle.fragment_numbers[specie_index]):
                optimized_coordinates[specie_index].append([])
                for i in range(len(particle.fragments[specie_index])):
                    optimized_coordinates[specie_index][fragment_number].append([int(file[coordinates_index+5+counter].split()[1]),float(file[coordinates_index+5+counter].split()[3]),float(file[coordinates_index+5+counter].split()[4]),float(file[coordinates_index+5+counter].split()[5])]) 
                    counter+=1
        for specie_index in range(particle.num_species):
            optimized_particle_coordinates.append([])
            for fragment_number in range(particle.fragment_numbers[specie_index]):
                optimized_particle_coordinates[specie_index].append([])
                optimized_particle_coordinates[specie_index][fragment_number].append(fragment_center_of_mass(optimized_coordinates[specie_index][fragment_number]))
                optimized_particle_coordinates[specie_index][fragment_number].append(angles_from_fragment(particle.fragments[specie_index],fragment_alignment(optimized_coordinates[specie_index][fragment_number])))
        
        return(optimized_particle_coordinates, energy)

def read_gaussian_file_with_align(particle):
    Normal_termination = True
    optimized_coordinates = []
    optimized_particle_coordinates = []
    with open(particle.name+'.log', 'r') as f:
        file = f.read()
        file = file.split('\n')
        if 'Error' in file[-5]:
            Normal_termination = False
        coordinates_index=len(file)-1-file[::-1].index('                          Input orientation:                          ')
        energy_index = -1
        for i in file[::-1]:
            energy_index+=1
            if 'SCF Done:' in i:
                energy_index=len(file)-1-energy_index
                break
        energy = float(file[energy_index].split()[4])
        counter=0
        for specie_index in range(particle.num_species):
            optimized_coordinates.append([])
            for fragment_number in range(particle.fragment_numbers[specie_index]):
                optimized_coordinates[specie_index].append([])
                for i in range(len(particle.fragments[specie_index])):
                    optimized_coordinates[specie_index][fragment_number].append([int(file[coordinates_index+5+counter].split()[1]),float(file[coordinates_index+5+counter].split()[3]),float(file[coordinates_index+5+counter].split()[4]),float(file[coordinates_index+5+counter].split()[5])]) 
                    counter+=1
        translation_shift = fragment_center_of_mass(optimized_coordinates[particle.Center_molecule[0]][particle.Center_molecule[1]])
        #trial_rotation=fragment_rotation_xyz(fragment_alignment(optimized_coordinates[particle.Center_molecule[0]][particle.Center_molecule[1]]),np.pi/4,np.pi/4,np.pi/4)
        global_rotation = angles_from_fragment(fragment_alignment(optimized_coordinates[particle.Center_molecule[0]][particle.Center_molecule[1]]),particle.fragments[particle.Center_molecule[0]])
        #global_rotation[0],global_rotation[1],global_rotation[2]=global_rotation[0]-np.pi/4,global_rotation[1]-np.pi/4,global_rotation[2]-np.pi/4
        #print(global_rotation)
        for specie_index in range(particle.num_species):
            for fragment_number in range(particle.fragment_numbers[specie_index]):
                    optimized_coordinates[specie_index][fragment_number]=fragment_shift(optimized_coordinates[specie_index][fragment_number], translation_shift[0], translation_shift[1], translation_shift[2])
                    optimized_coordinates[specie_index][fragment_number]=fragment_rotation_xyz(optimized_coordinates[specie_index][fragment_number], global_rotation[0], global_rotation[1], global_rotation[2])
        test = []                                       
        for i in optimized_coordinates:
            for j in i:
                for w in j:
                    test.append(w)
        #print(string_from_fragment(test))
        
        for specie_index in range(particle.num_species):
            optimized_particle_coordinates.append([])
            for fragment_number in range(particle.fragment_numbers[specie_index]):
                optimized_particle_coordinates[specie_index].append([])
                if [specie_index, fragment_number] == particle.Center_molecule:
                    optimized_particle_coordinates[specie_index][fragment_number].append([0.0,0.0,0.0])
                    optimized_particle_coordinates[specie_index][fragment_number].append([0.0,0.0,0.0])
                else:
                    optimized_particle_coordinates[specie_index][fragment_number].append(fragment_center_of_mass(optimized_coordinates[specie_index][fragment_number]))
                    optimized_particle_coordinates[specie_index][fragment_number].append(angles_from_fragment(particle.fragments[specie_index],fragment_alignment(optimized_coordinates[specie_index][fragment_number])))
        
        return(optimized_particle_coordinates, energy, Normal_termination)

    
def read_gaussian_file_with_align_no_rotation(particle):
    Normal_termination = True
    optimized_coordinates = []
    optimized_particle_coordinates = []
    with open(particle.name+'.log', 'r') as f:
        file = f.read()
        file = file.split('\n')
        if 'Error' in file[-5]:
            Normal_termination = False
        coordinates_index=len(file)-1-file[::-1].index('                          Input orientation:                          ')
        energy_index = -1
        for i in file[::-1]:
            energy_index+=1
            if 'SCF Done:' in i:
                energy_index=len(file)-1-energy_index
                break
        energy = float(file[energy_index].split()[4])
        counter=0
        for specie_index in range(particle.num_species):
            optimized_coordinates.append([])
            for fragment_number in range(particle.fragment_numbers[specie_index]):
                optimized_coordinates[specie_index].append([])
                for i in range(len(particle.fragments[specie_index])):
                    optimized_coordinates[specie_index][fragment_number].append([int(file[coordinates_index+5+counter].split()[1]),float(file[coordinates_index+5+counter].split()[3]),float(file[coordinates_index+5+counter].split()[4]),float(file[coordinates_index+5+counter].split()[5])]) 
                    counter+=1
        translation_shift = fragment_center_of_mass(optimized_coordinates[particle.Center_molecule[0]][particle.Center_molecule[1]])
        #trial_rotation=fragment_rotation_xyz(fragment_alignment(optimized_coordinates[particle.Center_molecule[0]][particle.Center_molecule[1]]),np.pi/4,np.pi/4,np.pi/4)
        global_rotation = angles_from_fragment(fragment_alignment(optimized_coordinates[particle.Center_molecule[0]][particle.Center_molecule[1]]),particle.fragments[particle.Center_molecule[0]])
        #global_rotation[0],global_rotation[1],global_rotation[2]=global_rotation[0]-np.pi/4,global_rotation[1]-np.pi/4,global_rotation[2]-np.pi/4
        #print(global_rotation)
        for specie_index in range(particle.num_species):
            for fragment_number in range(particle.fragment_numbers[specie_index]):
                    optimized_coordinates[specie_index][fragment_number]=fragment_shift(optimized_coordinates[specie_index][fragment_number], translation_shift[0], translation_shift[1], translation_shift[2])
                    optimized_coordinates[specie_index][fragment_number]=fragment_rotation_xyz(optimized_coordinates[specie_index][fragment_number], global_rotation[0], global_rotation[1], global_rotation[2])
        #test = []                                       
        #for i in optimized_coordinates:
        #    for j in i:
        #        for w in j:
        #            test.append(w)
        #print(string_from_fragment(test))
        full_optimized_fragment_coordinates =[]
        
        for specie_index in range(particle.num_species):
            optimized_particle_coordinates.append([])
            full_optimized_fragment_coordinates.append([])
            for fragment_number in range(particle.fragment_numbers[specie_index]):
                optimized_particle_coordinates[specie_index].append([])
                #full_optimized_fragment_coordinates[specie_index].append([])
                optimized_particle_coordinates[specie_index][fragment_number].append(fragment_center_of_mass(optimized_coordinates[specie_index][fragment_number]))
                full_optimized_fragment_coordinates[specie_index].append(fragment_alignment(optimized_coordinates[specie_index][fragment_number]))
        
        return(full_optimized_fragment_coordinates, optimized_particle_coordinates, energy, Normal_termination)
def grouping_and_updating(Swarm):
    energy=[]
    for i in Swarm:
        a, b, c, d=read_gaussian_file_with_align_no_rotation(i)
        i.fragment_coords = b
        i.all_fragments=a
        if i.Best_particle_results == None:
            i.Best_particle_results=[c, b]
        elif i.Best_particle_results[0] > c:
            i.Best_particle_results=[c, b]
        energy.append((c,d))
    energy=enumerate(energy)
    energy=list(energy)
    energy=sorted(energy, key = operator.itemgetter(1))
    shifted_energy=[((n[1][0]-energy[0][1][0])*627.5, n[0], n[1][1]) for n in energy]
    
    if Swarm[energy[0][0]].Best_swarm_results==None:
        for i in Swarm:
            i.Best_swarm_results=[energy[0][1][0], Swarm[energy[0][0]].fragment_coords]
    
    energy_previous=-1
    structures_of_one_type=[]
    for i in shifted_energy:
        if energy_previous - i[0] < -0.2:
            #print('New_structure_type: ', i[1], i[2])
            structures_of_one_type.append([i])
            energy_previous=i[0]
        else:
            structures_of_one_type[-1].append(i)
            
    corrected_structures_of_one_type=[]
    for same_structures in structures_of_one_type:
        proposed_clustering=[]
        for i in range(len(same_structures)):
            specie_other=[]
            specie_same=[]
            for j in range(len(same_structures)):
                if RMS_distance_evaluation(Swarm[same_structures[i][1]],Swarm[same_structures[j][1]])*other_RMS_distance_evaluation(Swarm[same_structures[i][1]],Swarm[same_structures[j][1]])>0.5:
                    specie_other.append(j)
                else: specie_same.append(j)
            proposed_clustering.append((specie_same, specie_other))
        new_groups=[]
        for i in proposed_clustering:
            cluster = True
            for same in i[0]:
                ans = True
                for k in i[0]:
                    ans*= k in proposed_clustering[k][0]
                if not bool(ans):
                    cluster=False
                    break
            if cluster:
                new_groups.append(i[0])
        newlist=[]
        for j in new_groups:
            if j not in newlist:
                newlist.append(j)
        #print(newlist)
        total_length = 0
        for i in newlist:
            total_length += len(i)
        if total_length == len(same_structures):
            for i in newlist:
                new_group=[]
                for j in i:
                    new_group.append(same_structures[j])
                corrected_structures_of_one_type.append(new_group)
        else: corrected_structures_of_one_type.append(same_structures)
    return(corrected_structures_of_one_type)
def coordinate_update(particle, w=0.7, c1=0.6, c2=0.7,  Discrete_randomization=True):
    if particle.Best_swarm_results == None:
        return
    else:
        for i in range(particle.num_species):
            for j in range(particle.fragment_numbers[i]):
                if Discrete_randomization:
                    rand1=np.random.choice([0,1])
                    rand2=np.random.choice([0,1])
                else:
                    rand1=np.random.random()
                    rand2=np.random.random()
                particle.velocity[i][j][0] = w*particle.velocity[i][j][0]+c1*rand1*(particle.Best_particle_results[1][i][j][0][0]-particle.fragment_coords[i][j][0][0])+c2*rand2*(particle.Best_swarm_results[1][i][j][0][0]-particle.fragment_coords[i][j][0][0])
                particle.velocity[i][j][1] = w*particle.velocity[i][j][1]+c1*rand1*(particle.Best_particle_results[1][i][j][0][1]-particle.fragment_coords[i][j][0][1])+c2*rand2*(particle.Best_swarm_results[1][i][j][0][1]-particle.fragment_coords[i][j][0][1])
                particle.velocity[i][j][2] = w*particle.velocity[i][j][2]+c1*rand1*(particle.Best_particle_results[1][i][j][0][2]-particle.fragment_coords[i][j][0][2])+c2*rand2*(particle.Best_swarm_results[1][i][j][0][2]-particle.fragment_coords[i][j][0][2])
                particle.fragment_coords[i][j][0][0] += particle.velocity[i][j][0]
                particle.fragment_coords[i][j][0][1] += particle.velocity[i][j][1]
                particle.fragment_coords[i][j][0][2] += particle.velocity[i][j][2]
        return
class Particle:
    def __init__(self, name, fragment_numbers, fragments, Box_size=5, Line=0.5, Plane=0.5, fragment_separation=0.4, Is_Line=False, Is_Plane=False, Is_Solvatation=False, Center_molecule=None):
        self.name = name
        self.Center_molecule=Center_molecule
        self.num_species=len(fragment_numbers)
        self.number_fragments=sum(fragment_numbers)
        self.fragment_numbers=fragment_numbers
        self.Box_size=Box_size
        self.Line=Line
        self.Plane=Plane
        self.fragment_separation=fragment_separation
        self.Is_Line = Is_Line
        self.Is_Plane = Is_Plane
        self.Is_Solvatation = Is_Solvatation
        self.fragments=[]
        self.fragments_radii=[]
        self.fragment_coords=[]
        self.Best_swarm_results=None
        self.Best_particle_results=None
        self.velocity=[]
        for i in range(self.num_species):
            self.velocity.append(np.zeros((self.fragment_numbers[i],3)))
                
        self.all_fragments=[]
        if self.Center_molecule == None:
            self.Center_molecule = [len(fragment_numbers)-1, fragment_numbers[-1]-1]
        for number_of_fragment in self.fragment_numbers:
            self.fragment_coords.append([])
            for i in range(number_of_fragment):
                self.fragment_coords[-1].append([np.zeros(3),np.zeros(3)])
        for i in fragments:
            self.fragments.append(i)
            self.fragments_radii.append(fragment_radii(i))
        for specie_index in range(len(self.fragment_numbers)):
            self.all_fragments.append([])
            for fragment_number in range(self.fragment_numbers[specie_index]):
                self.all_fragments[-1].append(self.fragments[specie_index])
    def check_overlap(self):
        overlaping_fragments=[]
        for index_of_specie in range(self.num_species):
            for number_of_fragment in range(self.fragment_numbers[index_of_specie]):
                restriction = self.fragments_radii[index_of_specie]*2+self.fragment_separation
                if number_of_fragment+1==self.fragment_numbers[index_of_specie]:
                    pass
                else:
                    for number_of_fragment_2 in range(number_of_fragment+1,self.fragment_numbers[index_of_specie]):
                        dist=np.sqrt(((self.fragment_coords[index_of_specie][number_of_fragment][0][0]-self.fragment_coords[index_of_specie][number_of_fragment_2][0][0])**2+
                              (self.fragment_coords[index_of_specie][number_of_fragment][0][1]-self.fragment_coords[index_of_specie][number_of_fragment_2][0][1])**2+
                              (self.fragment_coords[index_of_specie][number_of_fragment][0][2]-self.fragment_coords[index_of_specie][number_of_fragment_2][0][2])**2))
                        if dist<restriction:
                            overlaping_fragments.append(([index_of_specie,number_of_fragment],[index_of_specie,number_of_fragment_2])) #[([specie,fragment],[specie,fragment]),...]
            if index_of_specie+1==self.num_species:
                pass
            else:
                for number_of_fragment in range(self.fragment_numbers[index_of_specie]):
                    for index_of_specie_2 in range(index_of_specie+1, self.num_species):
                        restriction = self.fragments_radii[index_of_specie]+self.fragments_radii[index_of_specie_2]+self.fragment_separation
                        for number_of_fragment_2 in range(self.fragment_numbers[index_of_specie_2]):
                            #print("I'm checking:", ([index_of_specie,number_of_fragment],[index_of_specie_2,number_of_fragment_2]), 'with DR: ', restriction)
                            dist=np.sqrt(((self.fragment_coords[index_of_specie][number_of_fragment][0][0]-self.fragment_coords[index_of_specie_2][number_of_fragment_2][0][0])**2+
                                  (self.fragment_coords[index_of_specie][number_of_fragment][0][1]-self.fragment_coords[index_of_specie_2][number_of_fragment_2][0][1])**2+
                                  (self.fragment_coords[index_of_specie][number_of_fragment][0][2]-self.fragment_coords[index_of_specie_2][number_of_fragment_2][0][2])**2))
                            #print("Distanse is: ", dist)
                            if dist<restriction:
                                overlaping_fragments.append(([index_of_specie,number_of_fragment],[index_of_specie_2,number_of_fragment_2]))
        return(bool(len(overlaping_fragments)), overlaping_fragments)
    def fix_overlap(self):
	
        overlaped, overlaping_fragments = self.check_overlap()
        if overlaped:
            pair = overlaping_fragments[np.random.choice(len(overlaping_fragments))]
            permutation = np.random.permutation(2)
            fragment1 = pair[permutation[0]]
            fragment2 = pair[permutation[1]]
            self.fragment_coords[fragment1[0]][fragment1[1]][0][0] = self.fragment_coords[fragment1[0]][fragment1[1]][0][0] - (self.fragment_coords[fragment2[0]][fragment2[1]][0][0]-self.fragment_coords[fragment1[0]][fragment1[1]][0][0])*(np.random.random())
            self.fragment_coords[fragment1[0]][fragment1[1]][0][1] = self.fragment_coords[fragment1[0]][fragment1[1]][0][1] - (self.fragment_coords[fragment2[0]][fragment2[1]][0][1]-self.fragment_coords[fragment1[0]][fragment1[1]][0][1])*(np.random.random())
            self.fragment_coords[fragment1[0]][fragment1[1]][0][2] = self.fragment_coords[fragment1[0]][fragment1[1]][0][2] - (self.fragment_coords[fragment2[0]][fragment2[1]][0][2]-self.fragment_coords[fragment1[0]][fragment1[1]][0][2])*(np.random.random())
            self.fix_overlap()
        else:
            pass
            #print('Structure is perfect!')
    def random_initialization(self):
        if self.Is_Line:
            for index_of_specie in range(self.num_species):
                for number_of_fragment in range(self.fragment_numbers[index_of_specie]):
                    self.fragment_coords[index_of_specie][number_of_fragment][0]=([self.Box_size*np.random.random()-self.Box_size/2,self.Line*np.random.random()-self.Line/2,self.Line*np.random.random()-self.Line/2])
                    self.fragment_coords[index_of_specie][number_of_fragment][1]=(2*np.pi*np.random.random(3)-np.pi)
                    self.all_fragments[index_of_specie][number_of_fragment]=fragment_rotation_xyz(copy.deepcopy(self.all_fragments[index_of_specie][number_of_fragment]),*self.fragment_coords[index_of_specie][number_of_fragment][1])
        elif self.Is_Plane:
            for index_of_specie in range(self.num_species):
                for number_of_fragment in range(self.fragment_numbers[index_of_specie]):
                    self.fragment_coords[index_of_specie][number_of_fragment][0]=([self.Box_size*np.random.random()-self.Box_size/2,self.Box_size*np.random.random()-self.Box_size/2,self.Plane*np.random.random()-self.Plane/2])
                    self.fragment_coords[index_of_specie][number_of_fragment][1]=(2*np.pi*np.random.random(3)-np.pi)
                    self.all_fragments[index_of_specie][number_of_fragment]=fragment_rotation_xyz(copy.deepcopy(self.all_fragments[index_of_specie][number_of_fragment]),*self.fragment_coords[index_of_specie][number_of_fragment][1])
        elif self.Is_Solvatation:
            for index_of_specie in range(self.num_species):
                for number_of_fragment in range(self.fragment_numbers[index_of_specie]):
                    if index_of_specie == self.num_species-1:
                        pass
                    else:
                        self.fragment_coords[index_of_specie][number_of_fragment][0]=(self.Box_size*np.random.random(3)-self.Box_size/2)
                        self.fragment_coords[index_of_specie][number_of_fragment][1]=(2*np.pi*np.random.random(3)-np.pi)
                    self.all_fragments[index_of_specie][number_of_fragment]=fragment_rotation_xyz(copy.deepcopy(self.all_fragments[index_of_specie][number_of_fragment]),*self.fragment_coords[index_of_specie][number_of_fragment][1])
        else:
            for index_of_specie in range(self.num_species):
                for number_of_fragment in range(self.fragment_numbers[index_of_specie]):
                    self.fragment_coords[index_of_specie][number_of_fragment][0]=(self.Box_size*np.random.random(3)-self.Box_size/2)
                    self.fragment_coords[index_of_specie][number_of_fragment][1]=(2*np.pi*np.random.random(3)-np.pi)
                    self.all_fragments[index_of_specie][number_of_fragment]=fragment_rotation_xyz(copy.deepcopy(self.all_fragments[index_of_specie][number_of_fragment]),*self.fragment_coords[index_of_specie][number_of_fragment][1])
    
    def print_coordinates(self):
        coords_final=str()
        for index, coords in enumerate(self.fragment_coords):
            for d, i in enumerate(coords):
                new_fragment=fragment_placement(copy.deepcopy(self.all_fragments[index][d]), i[0][0],i[0][1],i[0][2])
                coords_final+=string_from_fragment(new_fragment)
        return(coords_final)
    def get_full_coordinates(self):
        coords_final=[]
        for index, coords in enumerate(self.fragment_coords):
            for d, i in enumerate(coords):
                new_fragment=fragment_placement(copy.deepcopy(self.all_fragments[index][d]), i[0][0],i[0][1],i[0][2])
                for k in new_fragment:
                    coords_final.append(k)
        return(coords_final)