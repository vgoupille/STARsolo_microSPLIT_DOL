#! /bin/bash
#SBATCH --job-name=create_env_STARsolo_Sequencing_Plateform
#SBATCH --mail-type=ALL
#SBATCH --mail-user=valentin.goupille@univ-rennes.fr
. /local/env/envconda.sh

conda create -p /home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/5_environnements/env_STARsolo python=3.6 -y

conda activate /home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/5_environnements/env_STARsolo

conda install -c bioconda star=2.7.9 cufflinks=2.2.1 -y


STAR --version
cufflinks --version



