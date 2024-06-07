#!/cm/shared/apps/python/3.7.x-anaconda/bin/python
#python use scrublet  ##jk use python/3.8.x-anaconda

## USAGE NOTES
# You can run it, after loading the appropriate module, with:
#
#   ./scrublet_data_process.py
# 
# In case you get an error related to ^R when trying to execute it as above,
# you may need to pass this file through dos2unix. 
# 
# In your environment, you may need to run `module load dos2unix/gcc/7.3.4`. 
# Afterwards, run `dos2unix scrublet_data_process.py` then execute as above. 
#

#tar xfz pbmc8k_filtered_gene_bc_matrices.tar.gz
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt 
import numpy as np
import os 
import pandas as pd 
import sys
import multiprocessing
import h5py
import readh5

start_dir = "./ao/data_pipeline1/"
data_path = "cellbender"


def calculate_doublets(cur_dir, d):
    print("[X] Exploring: " + cur_dir)

    dirs = os.listdir(cur_dir)
    for h5file in dirs:
        if h5file.endswith('.h5'):
            ## Filenames
            matrix_file = os.path.join(cur_dir, h5file)
            results_file = os.path.join(start_dir, d, h5file + '.scrublet_results.csv')
            results_figure_file = os.path.join(start_dir, d, h5file + '.scrublet_results.png')
            results_umap_file = os.path.join(start_dir, d, h5file + '.scrublet_umap_results.png')


            print("[{}] Building counts_matrix".format(h5file))
            #custom read function that needs to be loaded preceding running this function
            counts_matrix = readh5.anndata_from_h5(matrix_file, False)

            print('[{}] Counts matrix shape: {} rows, {} columns'.format(h5file, counts_matrix.shape[0], counts_matrix.shape[1]))
            print('[{}] Number of genes in gene list: {}'.format(h5file, counts_matrix.shape[1]))

            num_cells = counts_matrix.shape[0]
            expected_doublet_rate = num_cells/125000 # linear approximation to provided table
            print('[{}] Using an expected_doublet_rate of {}'.format(h5file, expected_doublet_rate))
            print('[{}] Scrubbing'.format(h5file))
            scrub = scr.Scrublet(counts_matrix.X, expected_doublet_rate=expected_doublet_rate)

            doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                                        min_cells=3, 
                                                                        min_gene_variability_pctl=85,
                                                                        n_prin_comps=35)

            print('[{}] Calculated doublet scores and predicted doublets, saving results to {}'.format(h5file, results_file))

            df = pd.DataFrame({
                'doublet_score': doublet_scores,
                'predicted_doublet': predicted_doublets
            })
            df.to_csv(results_file, index=False)

            print('[{}] Saving figure to {}'.format(h5file, results_figure_file))
            fig, axs = scrub.plot_histogram()
            fig.tight_layout(pad=2)
            fig.suptitle(d[:-2], fontsize=16)
            plt.savefig(results_figure_file)

            print('[{}] Running UMAP...'.format(h5file))
            scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))

            ###plot doublet on 3d projection
            print('[{}] Saving UMAP to {}'.format(h5file, results_umap_file))
            fig, axs = scrub.plot_embedding('UMAP', order_points=True)
            fig.tight_layout(pad=2)
            fig.suptitle(d[:-2], fontsize=16)
            plt.savefig(results_umap_file)

            print('------------------------------------')
            print('[{}] Done ^.^'.format(h5file))


program_name = """
 _____ _____ ______ _   _____________ ___________ 
/  ___/  __ \\| ___ \\ | | | ___ \\ ___ \\  ___| ___ \\
\\ `--.| /  \/| |_/ / | | | |_/ / |_/ / |__ | |_/ /
 `--. \\ |    |    /| | | | ___ \\ ___ \\  __||    / 
/\\__/ / \\__/\\| |\\ \\| |_| | |_/ / |_/ / |___| |\\ \\ 
\\____/ \\____/\\_| \\_|\\___/\\____/\\____/\\____/\\_| \\_|
                       - 16 -                             
                                                  """
print(program_name)

dirs = os.listdir(start_dir)
dirs.sort()

#print("[X] Exploring directories: " + str(dirs))
processes = []
print('------------------------------------')
for d in dirs:
    #if d == "CELLRANGER":
    #    break
    cur_dir = os.path.join(start_dir, d, data_path)
    if os.path.exists(cur_dir):
        p = multiprocessing.Process(target=calculate_doublets, args=(cur_dir, d))
        processes.append(p)

for p in processes:
    p.start()

print("[X] Created " + str(len(processes)) + " processes")
