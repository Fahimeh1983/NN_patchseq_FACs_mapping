for ri in {1..100}

do  

jobid=$ri
jobname="imb_10000_500_"
echo '#!/bin/bash'>subjob.bash
echo '#PBS -q celltypes'>>subjob.bash
echo '#PBS -N some_name'>>subjob.bash
echo '#PBS -m a'>>subjob.bash
echo '#PBS -r n'>>subjob.bash
echo '#PBS -l nodes=1:ppn=16'>>subjob.bash
echo '#PBS -l mem=16g,walltime=200:00:00'>>subjob.bash
echo '#PBS -o /allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Manuscript_patchseq_2019/Taxonomy/logs/'${jobname}${jobid//./-}'.out'>>subjob.bash
echo '#PBS -j oe'>>subjob.bash
echo 'source activate mypython3'>>subjob.bash
echo 'cd /allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Manuscript_patchseq_2019/Taxonomy'>>subjob.bash
echo 'python -m  Keras_imbalanced_trainer  --batch_size 500 --epochs 10000 --data_dir /allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Manuscript_patchseq_2019/  --output_dir new_keras_models/ --dropout 0.6 --n_features 4020 --n_hidden 10 --n_test_cells 1464 --n_celltypes 93 --run_iter '$ri>>subjob.bash
echo '...'

sleep 1

wait

qsub subjob.bash

echo 'Job: '$jobid' '

done


rm subjob.bash
