# run filter_matrix_allpops.py script
module load python
python ../../scripts/filter_matrix_allpops.py

# run plot_r_allpops.py script
python ../../scripts/plot_r_allpops.py

# run plot_corrprob_allpops.py script
python ../../scripts/plot_corrprob_allpops.py

# run plot_joint_corr.py script
python ../../scripts/plot_joint_corr.py

# run run_annovar.sh script
bash ../../scripts/run_annovar.sh

# run combine_annovar.sh script
bash ../../scripts/combine_annovar.sh

# run plot_annotation_table.py script
python ../../scripts/plot_annotation_table.py

# run plot_annotation_contig.py script for each contiguous annotation set
python ../../scripts/plot_corr_annotations.py

# move all output plots to a new folder 
mkdir -p output_plots && mv *.pdf output_plots/
