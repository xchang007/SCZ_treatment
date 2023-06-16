# 1. Set up

# https://fmriprep.org/en/stable/docker.html
docker pull nipreps/fmriprep:latest


# 2. (Optional): Install fmriprep-docker wrapper 
conda activate PreProcess
python -m pip install --user --upgrade fmriprep-docker


# 3. S3_Update_PBS_JOB.m： To run fmriprep, create PBS_Job*.sh for each subject 


# 4. Create x_qsub.sh 


# 5. run x_qsub.sh 
sh x_qsub.sh /share/inspurStorage/home1/changxiao/Project/2021_SCZ_Cui/Code/PBS_Jobs_*



