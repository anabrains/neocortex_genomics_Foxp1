module load meme/5.3.3 
## cd ./sc-atac/aoa/gene_lists/fasta/fimo/
# Specify  directory w/ folders containing input files
input_dir="./sc-atac/aoa/gene_lists/fasta/fasta_seq"
save_dir="./sc-atac/aoa/gene_lists/fasta/fimo"

# motif file #downloaded from JASPAR
motif_file="./sc-atac/aoa/gene_lists/fasta/FOX_motifs_meme.txt"

# Loop through folders specifying overlapping degs & dars
for folder in "$input_dir"/*/; do
    ## get folder name
    folder_name=$(basename "$folder")
        echo "created the folder $folder_name"
        echo "also this was file path $folder"


    ## verify that targets.fa exists #aka genes in fa format
    if [ -e "$folder/target.fa" ]; then
        # Create output directory based on the current folder
        output_dir="$save_dir/${folder_name}_FIMO_output"
        mkdir -p "$output_dir"

        echo "ran mkdir called $output_dir"

        ## FIMO command
        fimo --oc "$output_dir" "$motif_file" "$folder"/target.fa
            ##    fimo --oc "$output_dir" --motif "$motif_file" "$folder"target.fa

    else
        echo "target.fa not found in $folder"
        ls "$folder" # for status
    fi
done