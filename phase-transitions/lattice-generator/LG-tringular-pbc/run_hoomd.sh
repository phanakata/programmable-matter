srun -p gpu --gres=gpu:1 --mem 1000 -n 4 --pty -t 600 /bin/bash 
singularity shell --nv software.simg

singularity exec software.simg python3 --version

srun --pty -p gpu -t 0-06:00 --mem 8000 --gres=gpu:1 /bin/bash 

srun -p test --pty --mem 500 -t 0-06:00 /bin/bash
