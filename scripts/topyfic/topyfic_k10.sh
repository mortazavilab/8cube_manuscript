#!/bin/sh
#SBATCH -A seyedam_lab
#SBATCH --cpus-per-task 1
#SBATCH --output=LDA_k10.out
#SBATCH --error=LDA_k10.err
#SBATCH --time=1:00:00
#SBATCH -J topyfic_control
#SBATCH --partition=standard

source ~/miniconda3/bin/activate scanpy_env

data=/share/crsp/lab/seyedam/share/igvf_pipeline/topyfic/Satellite_regulatory_only_normalized.h5ad
k=10
tissue=Satellite

for i in {0..99}
do
    scriptName=run_${i}
    curr=${scriptName}.sh
    echo '#!/bin/bash' > ${curr}
    echo '#SBATCH -A seyedam_lab' >> ${curr}
    echo '#SBATCH --cpus-per-task 10' >> ${curr}
    echo '#SBATCH --output=LDA-%J.out' >> ${curr}
    echo '#SBATCH --error=LDA-%J.err' >> ${curr}
    echo '#SBATCH --time=02:00:00' >> ${curr}
    echo '#SBATCH -J topyfic-%J' >> ${curr}
    echo '#SBATCH --partition=standard' >> ${curr}
    
    echo "source ~/miniconda3/bin/activate scanpy_env"
   
    echo "python3 topyfic_train.py -d ${data} -t ${tissue} -k ${k} -r ${i}" >> ${curr}
    
    chmod +x ${curr}
    sbatch ${scriptName}.sh
    
done