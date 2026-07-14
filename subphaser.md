# Install SubPhaser

Follow the instruction in [this](https://github.com/zhangrengang/SubPhaser) repository.
```sh
git clone https://github.com/zhangrengang/SubPhaser
cd SubPhaser

# Install
conda env create -f SubPhaser.yaml
conda activate SubPhaser
python setup.py install
```

# Run SubPhaser
Bash script need to run the SubPhaser analysis.
```sh
#!/bin/bash
#SBATCH --job-name=SubPhaser
#SBATCH -o ./logs/subphaser.%j.out
#SBATCH -e ./logs/subphaser.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=500:00:00
#SBATCH --mem-per-cpu=3800
#SBATCH --nice=1000
#SBATCH --partition=main

conda activate SubPhaser

subphaser -i out_JBAT_review4.FINAL_22scaf_rename.fa.gz -c sg.config -p 32 -max_memory 7.8G

conda deactivate
```
