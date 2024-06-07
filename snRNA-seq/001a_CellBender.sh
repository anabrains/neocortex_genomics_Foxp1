# load env
source activate CellBender

##should be able to run the CellBender now
cellbender remove-background \
                 --input raw_feature_bc_matrix.h5 \
                 --output ./ao/eao/cellBender_process/output/AO1_output.h5 \
                 --expected-cells 10000 \
                 --total-droplets-included 25000 \
                 --fpr 0.01 \
                 --epochs 150 --cuda
