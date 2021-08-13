import os
from subprocess import check_output
import time
import pickle

#reading variables from the Input.txt file

with open('Input.txt', 'r') as f:
    file = f.read()
file = file.split('\n')
file = [i.split() for i in file]

main_path = file[2][1]
script_file_header = file[3][1]
num_generations = int(file[4][1])
Generation_size = int(file[8][1])

if not os.path.exists("Generation.data"):

    os.chdir(main_path)
    os.system('python Swarm_creation.py')
    os.chdir(main_path+'/Generation_1_P3')
    job_list=[]
    for i in range(Generation_size):
        s=check_output('sbatch '+script_file_header+'_particle_'+str(i), shell=True)
        job_list.append(int(s.split()[3]))
    
    print("Job_list: ", job_list)

    os.chdir(main_path)
    time.sleep(300)
    jobs_run=2
    while jobs_run>0:
        jobs_run=0
        for i in job_list:
            try:
                s=check_output('squeue -j '+str(i), shell=True)
            except:
                print('Exception XXX')
                s='0'
            if 'R    ' in str(s) or 'PD ' in str(s) or 'CG    ' in str(s) or 'CF     ' in str(s):
                jobs_run+=1
        print(str(jobs_run)+' jobs are still running')
        time.sleep(300)

for i in range(num_generations-1):
    with open('Generation.data', 'rb') as f:
        Generation=pickle.load(f)
    if not os.path.exists(main_path+'/Generation_'+str(Generation)):
        os.system('python Swarm_presice_setup.py')
        os.chdir(main_path+'/Generation_'+str(Generation))
        job_list=[]
        for i in range(Generation_size):
            s=check_output('sbatch '+script_file_header+'_particle_'+str(i), shell=True)
            job_list.append(int(s.split()[3]))
    
        print("Job_list: ", job_list)

        os.chdir(main_path)
        time.sleep(300)
        jobs_run=2
        while jobs_run>0:
            jobs_run=0
            for i in job_list:
                try:
                    s=check_output('squeue -j '+str(i), shell=True)
                except:
                    print('Exception XXX')
                    s='0'
                if 'R    ' in str(s) or 'PD ' in str(s) or 'CG    ' in str(s) or 'CF     ' in str(s):
                    jobs_run+=1
            print(str(jobs_run)+' jobs are still running')
            time.sleep(300)


    os.system('python Swarm_raw_setup.py')
    with open('Generation.data', 'rb') as f:
        Generation=pickle.load(f)
    os.chdir(main_path+'/Generation_'+str(Generation)+'_P3')
    job_list=[]
    for i in range(Generation_size):
        s=check_output('sbatch '+script_file_header+'_particle_'+str(i), shell=True)
        job_list.append(int(s.split()[3]))
    
    print("Job_list: ", job_list)

    os.chdir(main_path)
    time.sleep(300)
    jobs_run=2
    while jobs_run>0:
        jobs_run=0
        for i in job_list:
            try:
                s=check_output('squeue -j '+str(i), shell=True)
            except:
                print('Exception XXX')
                s='0'
            if 'R    ' in str(s) or 'PD ' in str(s) or 'CG    ' in str(s) or 'CF     ' in str(s):
                jobs_run+=1
        print(str(jobs_run)+' jobs are still running')
        time.sleep(300)


    
    
