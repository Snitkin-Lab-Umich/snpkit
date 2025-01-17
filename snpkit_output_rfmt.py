import os
import shutil
import sys
import glob

"""
Reformatting output results from snpkit (python version)——creating the new folder structure for snpkit
"""

# Define a global variable for the error log
error_log = None

def check_required_dirs_and_files(base_path):
    """
    Check if the required directories and files exist i.e. snpkit ran successfully
    and that they are not empty.
    """
    core_results_pattern = '*_core_results'
    
    # Check for core results directory with unknown prefix
    core_results_dirs = glob.glob(os.path.join(base_path, core_results_pattern))
    if not core_results_dirs:
        print(f"\nDirectory ending in core_results not found in {base_path}. Did you finish running snpkit?")
        sys.exit(1)
    
    core_results_dir = core_results_dirs[0]  # Taking the first matching directory
    if not os.listdir(core_results_dir):
        print(f"\n{core_results_dir} is empty. Did you finish running snpkit?")
        sys.exit(1)
    else:
        print(f"Directory {core_results_dir} exists and is not empty.")
    
    # Define and check required subdirectories within the core results directory
    needed_dirs = ['data_matrix', 'gubbins']
    for dir_name in needed_dirs:
        dir_path = os.path.join(core_results_dir, dir_name)
        if not os.path.exists(dir_path):
            print(f"\n{dir_name} does not exist in {dir_path}. Did you finish running snpkit?")
            sys.exit(1)
        elif not os.listdir(dir_path):  # Check if the directory is empty
            print(f"\n{dir_name} is empty in {dir_path}. Did you finish running snpkit?")
            sys.exit(1)
        else:
            print(f"\n{dir_path} exists and is not empty.")
    return core_results_dir
    

def create_folders(base_path):
    """
    Creating all the folders needed first before files are moved from the results folder
    """
    snpkit_results_path = os.path.dirname(base_path)
    ref_genome_folder = os.path.join(snpkit_results_path,"snpkit_rfmt_results","ref_genome_filtered_positions")
    final_vcfs_folder = os.path.join(snpkit_results_path,"snpkit_rfmt_results","final_vcfs")
    depth_coverage_folder = os.path.join(snpkit_results_path,"snpkit_rfmt_results","depth_of_coverage_stats")
    variant_matrices_folder = os.path.join(snpkit_results_path,"snpkit_rfmt_results","variant_matrices")
    temp_variant_matrices_folder = os.path.join(snpkit_results_path,"snpkit_rfmt_results","variant_matrices", "temp_files")
    annotated_files_folder = os.path.join(snpkit_results_path, "snpkit_rfmt_results", "annotated_files")
    phylo_analysis_folder = os.path.join(snpkit_results_path,"snpkit_rfmt_results","phylo_analysis")
    gubbins_folder = os.path.join(snpkit_results_path,"snpkit_rfmt_results","phylo_analysis", "gubbins")
    genome_alignment_folder = os.path.join(snpkit_results_path,"snpkit_rfmt_results","phylo_analysis", "genome_alignment")


    try:
        os.makedirs(ref_genome_folder, exist_ok=True)
        os.makedirs(final_vcfs_folder, exist_ok=True)
        os.makedirs(depth_coverage_folder, exist_ok=True)
        os.makedirs(variant_matrices_folder, exist_ok=True)
        os.makedirs(temp_variant_matrices_folder, exist_ok=True)
        os.makedirs(annotated_files_folder, exist_ok=True)
        os.makedirs(phylo_analysis_folder, exist_ok=True)
        os.makedirs(gubbins_folder, exist_ok=True)
        os.makedirs(genome_alignment_folder, exist_ok=True)

        print(f"\nVariant calling folders created successfully at {snpkit_results_path}.")

        global error_log  # Use the global variable to set the error log path
        snpkit_rfmt_results = os.path.join(snpkit_results_path, "snpkit_rfmt_results")
        #os.makedirs(snpkit_rfmt_results, exist_ok=True)

        # Create the error log file in the same directory as snpkit_rfmt_results
        error_log = os.path.join(snpkit_rfmt_results, "error_log.txt")
        with open(error_log, "w") as f:
            f.write("Error Log Initialized\n")

    except Exception as e:
        print(f"Error creating folders: {str(e)}")
    


# def move_files(src, dst, files):
#     """
#     Move a list of files from source to destination
#     """
    
#     for file_name in files:
#         src_files = os.path.join(src, file_name)
#         dst_file = os.path.join(dst, file_name)
#         #print(f'Moving {src_file} to {dst_file}')
#         for src_file in glob.glob(src_files):
#             if os.path.exists(src_file):
#                 os.makedirs(dst, exist_ok=True)
#                 shutil.move(src_file, dst_file)
#                 print("Successfully moved {src_file} to {dst_file}")
#             else:
#                 print("Unable to move files: {src_file}")

# def move_files(src, dst, files):
#     """
#     Move a list of files (or wildcard patterns) from source to destination.
#     """
#     for file_name in files:
#         src_files = os.path.join(src, file_name)
#         matched_files = glob.glob(src_files)

#         if not matched_files:
#             #print(f"No files matched the pattern: {src_files}")
#             raise FileNotFoundError(f"\nNo files matched the pattern: {src_files}")

#         for src_file in matched_files:
#             if os.path.exists(src_file):
#                 os.makedirs(dst, exist_ok=True)
#                 dst_file = os.path.join(dst, os.path.basename(src_file))
#                 shutil.move(src_file, dst_file)
#                 print(f"\nSuccessfully moved {src_file} to {dst_file}")
#             else:
#                 src_file
#                 print(f"\nFile {src_file} does not exist. Skipping.")

def move_files(src, dst, files):
    """
    Move a list of files (or wildcard patterns) from source to destination.
    """
    global error_log  # Access the global error log variable

    for file_name in files:
        src_files = os.path.join(src, file_name)
        matched_files = glob.glob(src_files)

        if not matched_files:
            with open(error_log, "a") as error_file:
                print(f"\nNo files matched the pattern: {src_files}")
                error_file.write(f"\nNo files matched the pattern: {src_files} \n")
            continue

        for src_file in matched_files:
            if os.path.exists(src_file):
                os.makedirs(dst, exist_ok=True)
                dst_file = os.path.join(dst, os.path.basename(src_file))
                shutil.move(src_file, dst_file)
                print(f"\nSuccessfully moved {src_file} to {dst_file}")
            else:
                with open(error_log, "a") as error_file:
                    print(f"\nFile {src_file} does not exist. Skipping.")
                    error_file.write(f"\nFile {src_file} does not exist. Skipping.\n")

            
def organize_data_matrix(base_path, core_results_path):
    """
    Organize files in the data_matrix folder
    """
    src_path = os.path.join(base_path,  core_results_path, 'data_matrix')
    snpkit_results_path = os.path.dirname(base_path)

    # Check that required subdirectories exist and have files
    # required_subdirs = ['Functional_annotation_results', 'logs', 'matrices', 'plots', 'snpEff_results']
    # for subdir in required_subdirs:
    #     subdir_path = os.path.join(src_path, subdir)
    #     if not os.path.exists(subdir_path) or not os.listdir(subdir_path):
    #         print(f"Required subdirectory {subdir_path} does not exist or is empty.")
    #         sys.exit(1)

    # Move files from core_results/data_matrix/Functional_annotation_results to ref_genome_filtered_positions/
    dst_structure_func_ann = {
        os.path.join(snpkit_results_path, "snpkit_rfmt_results", 'ref_genome_filtered_positions'): [
            "Functional_class_filter_positions.txt",
            "inexact_repeat_region_positions.txt",
            "phage_region_positions.txt",
            "repeat_region_positions.txt"
        ]
        # ,
        # os.path.join(base_path, 'ref_genome_filtered_positions'): [
        #     "phage_region_positions.txt"
        # ],
        # os.path.join(base_path, 'ref_genome_filtered_positions'): [
        #     "repeat_region_positions.txt"
        # ]
    }
    for dst, files in dst_structure_func_ann.items():
        move_files(os.path.join(src_path, 'Functional_annotation_results'), dst, files)

    # Move logs directory to the same level as core_snp_consensus
    # logs_src = os.path.join(src_path, 'logs')
    # logs_dst = os.path.join(base_path, 'logs')
    # if os.path.exists(logs_src):
    #     shutil.move(logs_src, logs_dst)

    # Move files from core_results/data_matrix/matrices/ to variant_matrices
    dst_structure_matrices = {
        os.path.join(snpkit_results_path, "snpkit_rfmt_results",'variant_matrices',  'temp_files'): [
            "All_indel_label_final_ordered_sorted.txt",
            "All_label_final_ordered_sorted.txt",
            "unique_indel_positions_file",
            "unique_positions_file"
        ],
        os.path.join(snpkit_results_path, "snpkit_rfmt_results", 'variant_matrices'): [
            "Indel_matrix_allele.tsv",
            "Indel_matrix_code.tsv",
            "Indel_matrix_code_unmasked.tsv",
            "SNP_matrix_allele_new.tsv",
            "SNP_matrix_allele_unmasked.tsv",
            "SNP_matrix_code.tsv",
            "SNP_matrix_code_unmasked.tsv"
        ]
    }
    for dst, files in dst_structure_matrices.items():
        move_files(os.path.join(src_path, 'matrices'), dst, files)

def organize_gubbins(base_path, core_results_path):
    """
    Organize files in the gubbins folder.
    """
    src_path = os.path.join(base_path, core_results_path, 'gubbins')
    snpkit_results_path = os.path.dirname(base_path)

    dst_structure_parent_dir = {
        os.path.join(snpkit_results_path, "snpkit_rfmt_results",'phylo_analysis', 'genome_alignment'): [
            "*_allele_unmapped_gubbins_masked.fa",
            "*_allele_unmapped.fa"
        ],
        os.path.join(snpkit_results_path, "snpkit_rfmt_results", 'phylo_analysis', 'gubbins'): [
            "*_allele_unmapped.filtered_polymorphic_sites.fasta",
            "*_allele_unmapped.filtered_polymorphic_sites.phylip",
            "*_allele_unmapped.final_tree.tre",
            "*_allele_unmapped_masked_recomb_positions.txt",
            "*_allele_unmapped.node_labelled.final_tree.tre",
            "*_allele_unmapped.per_branch_statistics.csv",
            "*_allele_unmapped.recombination_predictions.embl",
            "*_allele_unmapped.recombination_predictions.gff",
            "*_allele_unmapped.summary_of_snp_distribution.vcf"
        ]
    }

    for dst, files in dst_structure_parent_dir.items():
        move_files(src_path, dst, files)
    
    src_path_iqtree = os.path.join(base_path, core_results_path, 'gubbins', 'iqtree_masked_wga')

    dst_structure_iqtree_masked_wga = {
        os.path.join(snpkit_results_path, "snpkit_rfmt_results", 'phylo_analysis', 'tree'): [
            "*_allele_unmapped.bionj",
            "*_allele_unmapped.ckp.gz",
            "*_allele_unmapped.contree",
            "*_allele_unmapped.iqtree",
            "*_allele_unmapped.log",
            "*_allele_unmapped.mldist",
            "*_allele_unmapped.parstree",
            "*_allele_unmapped.splits.nex",
            "*_allele_unmapped.treefile"
        ]
    }
    
    for dst, files in dst_structure_iqtree_masked_wga.items():
        move_files(src_path_iqtree, dst, files)

def get_samples(base_path, core_results_path):
    """
    Get list of sample names from the snpkit_results directory.

    Args:
        base_path (str): Path to the snpkit_results directory.
        core_results_path (str): Path to the *_core_results folder.

    Returns:
        list: List of sample names.
    """
    # Define folders to exclude
    excluded_folders = {
        os.path.basename(core_results_path),
        "core_temp_dir",
        "Logs",
        "temp_jobs",
        "snpkit_rfmt_results",
    }

    sample_names = [] # Initialize the list of samples

    # Iterate through the contents of the base path
    for sample_folder in os.listdir(base_path):
        sample_folder_path = os.path.join(base_path, sample_folder)

        # Check if the item is a directory and not in the excluded list
        if os.path.isdir(sample_folder_path) and sample_folder not in excluded_folders:
            # Define paths to the subdirectories for this sample
            sample_stats_results = os.path.join(sample_folder_path, f"{sample_folder}_stats_results")
            sample_vcf_results = os.path.join(sample_folder_path, f"{sample_folder}_vcf_results")

            # Check if the required paths exist
            if os.path.isdir(sample_stats_results) and os.path.isdir(sample_vcf_results):
                # Append the sample to the list if conditions are met
                sample_names.append(sample_folder)

                # Create the necessary folders for this sample
                snpkit_rfmt_results_path = os.path.dirname(base_path)
                os.makedirs(os.path.join(snpkit_rfmt_results_path, "snpkit_rfmt_results", "annotated_files", sample_folder), exist_ok=True)
                os.makedirs(os.path.join(snpkit_rfmt_results_path, "snpkit_rfmt_results", "depth_of_coverage_stats", sample_folder), exist_ok=True)
                os.makedirs(os.path.join(snpkit_rfmt_results_path, "snpkit_rfmt_results", "final_vcfs", sample_folder), exist_ok=True)

    print(f"\nSamples in this snpkit run are: {sample_names}")
    return sample_names


def organize_samples_folder(samples_list, base_path):
    """
    Organize files in the samples folder
    """
    snpkit_results_path = os.path.dirname(base_path)
    
    for sample in samples_list:
        src_path_depth_cov = os.path.join(base_path, sample, f"{sample}_stats_results")

        # Move depth of coverage stats files
        dst_structure_depth_cov = {
            os.path.join(snpkit_results_path, "snpkit_rfmt_results", "depth_of_coverage_stats", f"{sample}"): [
                f"{sample}_alignment_stats",
                f"{sample}_depth_of_coverage",
                f"{sample}_depth_of_coverage.sample_summary"
            ],
        }

        for dst, files in dst_structure_depth_cov.items():
            move_files(src_path_depth_cov, dst, files)
    
        # Move sample specific  results i.e. raw, final indel and snp 
        src_path_vcf_results = os.path.join(base_path, sample, f"{sample}_vcf_results")

        dst_structure_vcf_results = {
            os.path.join(snpkit_results_path, "snpkit_rfmt_results", "final_vcfs", f"{sample}"): [
                f"{sample}_aln_mpileup_raw.vcf", 
                f"{sample}_aln_mpileup_raw.vcf_5bp_indel_removed.vcf.gz",
                f"{sample}_aln_mpileup_raw.vcf_5bp_indel_removed.vcf.gz.tbi",
                f"{sample}_aln_mpileup_raw.vcf_indel.vcf",                     
                f"{sample}_aln_mpileup_raw.vcf_indel.vcf.gz",                 
                f"{sample}_aln_mpileup_raw.vcf_indel.vcf.gz.tbi",              
                f"{sample}_aln_mpileup_raw.vcf_indel.vcf.idx",                
                f"{sample}_aln_mpileup_raw.vcf.gz",
                f"{sample}_aln_mpileup_raw.vcf.gz.tbi", 
                f"{sample}_aln_mpileup_raw.vcf.idx",
                f"{sample}_filter2_final.vcf_no_proximate_snp.vcf", 
                f"{sample}_filter2_final.vcf_no_proximate_snp.vcf.gz", 
                f"{sample}_filter2_final.vcf_no_proximate_snp.vcf.gz.tbi",
                f"{sample}_filter2_indel_final.vcf",
                f"{sample}_filter2_indel_final.vcf.gz", 
                f"{sample}_filter2_indel_final.vcf.gz.tbi"
            ],
        }

        for dst, files in dst_structure_vcf_results.items():
            move_files(src_path_vcf_results, dst, files)
        
        # Move annotated files of raw, snp and indel files 
        dst_structure_anno_vcf_results = {
            os.path.join(snpkit_results_path, "snpkit_rfmt_results", "annotated_files", f"{sample}"): [
                f"{sample}_aln_mpileup_raw.vcf_ANN.csv", 
                f"{sample}_aln_mpileup_raw.vcf_ANN.genes.txt",
                f"{sample}_aln_mpileup_raw.vcf_ANN.vcf", 
                f"{sample}_aln_mpileup_raw.vcf_ANN.vcf.gz",
                f"{sample}_aln_mpileup_raw.vcf_ANN.vcf.gz.tbi", 
                f"{sample}_filter2_final.vcf_no_proximate_snp.vcf_ANN.csv",
                f"{sample}_filter2_final.vcf_no_proximate_snp.vcf_ANN.genes.txt",
                f"{sample}_filter2_final.vcf_no_proximate_snp.vcf_ANN.vcf",
                f"{sample}_filter2_final.vcf_no_proximate_snp.vcf_ANN.vcf.gz",
                f"{sample}_filter2_final.vcf_no_proximate_snp.vcf_ANN.vcf.gz.tbi",
                f"{sample}_filter2_indel_final.vcf_ANN.csv",
                f"{sample}_filter2_indel_final.vcf_ANN.genes.txt",
                f"{sample}_filter2_indel_final.vcf_ANN.vcf",
                f"{sample}_filter2_indel_final.vcf_ANN.vcf.gz",
                f"{sample}_filter2_indel_final.vcf_ANN.vcf.gz.tbi"
            ],
        }

        for dst, files in dst_structure_anno_vcf_results.items():
            move_files(src_path_vcf_results, dst, files)

        # sample_vcf_results_path = os.path.join(base_path, sample, '{sample}_vcf_results')
        # raw_files =  ["{sample}_aln_mpileup_raw.vcf", "{sample}_aln_mpileup_raw.vcf.gz","{sample}_aln_mpileup_raw.vcf.gz.tbi", "{sample}_aln_mpileup_raw.vcf.idx"]
        # snp_files = ["{sample}_filter2_final.vcf_no_proximate_snp.vcf", "{sample}_filter2_final.vcf_no_proximate_snp.vcf.gz", "{sample}_filter2_final.vcf_no_proximate_snp.vcf.gz.tbi"]
        # indel_files = ["{sample}_filter2_indel_final.vcf","{sample}_filter2_indel_final.vcf.gz", "{sample}_filter2_indel_final.vcf.gz.tbi"]
        

        #sample_coverage_path = os.path.join(base_path, sample, '{sample}_stats_results')
        #coverage_files = ["{sample}_alignment_stats", "{sample}_depth_of_coverage", "{sample}_depth_of_coverage.sample_sumamry"]
        #move_files(sample_coverage_path,depth_coverage_folder, files)

def rename_files(samples_list, base_path):
    """
    Rename files in the final_vcfs and annotated_files folders based on specific conditions.
    """
    snpkit_results_path = os.path.join(os.path.dirname(base_path), "snpkit_rfmt_results")
    folders_to_rename = ["final_vcfs", "annotated_files"]

    for folder in folders_to_rename:
        folder_path = os.path.join(snpkit_results_path, folder)
        
        for sample in samples_list:
            sample_folder_path = os.path.join(folder_path, sample)  # Path to the sample's folder
            
            if not os.path.exists(sample_folder_path):
                raise FileNotFoundError(f"Sample folder {sample_folder_path} does not exist.")
            
            # Iterate over the files in the sample's folder
            for file_name in os.listdir(sample_folder_path):
                file_path = os.path.join(sample_folder_path, file_name)
                
                # Ensure it's a file, not a directory
                if os.path.isfile(file_path):
                    # Apply renaming rules
                    if file_name.startswith(f"{sample}_aln_mpileup_raw.vcf_indel"):
                        new_name = file_name.replace(f"{sample}_aln_mpileup_raw.vcf_indel", f"{sample}_raw_indel", 1)
                    elif file_name.startswith(f"{sample}_aln_mpileup_raw.vcf_5bp_indel_removed"):
                        new_name = file_name.replace(f"{sample}_aln_mpileup_raw.vcf_5bp_indel_removed", f"{sample}_raw_snp_5bp_indel_removed", 1)
                    elif file_name.startswith(f"{sample}_aln_mpileup_raw"):
                        new_name = file_name.replace(f"{sample}_aln_mpileup_raw", f"{sample}_raw_snp", 1)
                    elif file_name.startswith(f"{sample}_filter2_final.vcf_no_proximate_snp"):
                        new_name = file_name.replace(f"{sample}_filter2_final.vcf_no_proximate_snp", f"{sample}_final_snp", 1)
                    elif file_name.startswith(f"{sample}_filter2_indel_final"):
                        new_name = file_name.replace(f"{sample}_filter2_indel_final", f"{sample}_final_indel", 1)
                    else:
                        continue  # Skip files that don't match any rule

                    # Full path for the new file name
                    new_file_path = os.path.join(sample_folder_path, new_name)
                    
                    # Rename the file
                    os.rename(file_path, new_file_path)
                    print(f"\nRenamed {file_path} to {new_file_path}")

def files_moved():
    global error_log
    with open(error_log, 'r') as f:
        lines = f.readlines()
        
    if len(lines) > 1:
        print(f"\nSome files are missing and/or not moved to the new spkit output results directory. More info can be found here: {error_log}.")
    else:
        print("\nDirectory successfully reformatted!")
        os.remove(error_log)

def main(base_path):
    core_results_path = check_required_dirs_and_files(base_path)
    create_folders(base_path)
    organize_data_matrix(base_path, core_results_path)
    organize_gubbins(base_path, core_results_path)
    samples_list = get_samples(base_path, core_results_path)
    organize_samples_folder(samples_list, base_path)
    rename_files(samples_list, base_path)
    files_moved()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python snpkit_output_file_rfmt.py <PathToSnpkitResults>") # Fix the usage command
        sys.exit(1)

    usr_base_path = sys.argv[1]
    base_path = os.path.normpath(usr_base_path)
    main(base_path)