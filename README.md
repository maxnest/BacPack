# BacPack

## Step 1: Creating a config file inside a docker container

Usage: ./ConfigGen.sh --config <string> --species_tag <string> [--fasta <string>] [--r1_fq <string>] [--r2_fq <string>] [--long_fq <string>] [--medaka_model <string>] --output_dir <string> --ref_nucl_fasta_path <string> --min_len <string> 

	--config	The name of the configuration file in which the parameters for launching the pipeline will be written. As an example, 'BacPack_config.yaml'
	--species_tag	Species and/or strain name to be used when generating analysis results. As an example, 'Bt_st1402'
	--fasta Full path to the fasta file with assembled genome
	--r1_fq	Full path to the archive with raw forward (R1) short paired-end reads
	--r2_fq Full path to the archive with raw reverse (R2) short paired-end reads 
	--output_dir The name of the directory where all the results of the pipeline will be written
	--ref_nucl_fasta_path	Full path to the directory containing fasta-files with genome assemblies that will be considered as references
	--min_len	Minimum length of assembled sequences
	--long_fq	Full path to the archive with raw long reads
	--medaka_model	The name of the model that will be used to create consensus based on the results of the assembly of long reads. Available options: r103_fast_g507 r103_hac_g507 r103_min_high_g345 r103_min_high_g360 r103_prom_high_g360 r103_sup_g507 r1041_e82_260bps_fast_g632 r1041_e82_260bps_hac_g632 r1041_e82_260bps_sup_g632 r1041_e82_400bps_fast_g615 r1041_e82_400bps_fast_g632 r1041_e82_400bps_hac_g615 r1041_e82_400bps_hac_g632 r1041_e82_400bps_sup_g615 r104_e81_fast_g5015 r104_e81_hac_g5015 r104_e81_sup_g5015 r104_e81_sup_g610 r10_min_high_g303 r10_min_high_g340 r941_e81_fast_g514 r941_e81_hac_g514 r941_e81_sup_g514 r941_min_fast_g303 r941_min_fast_g507 r941_min_hac_g507 r941_min_high_g303 r941_min_high_g330 r941_min_high_g340_rle r941_min_high_g344 r941_min_high_g351 r941_min_high_g360 r941_min_sup_g507 r941_prom_fast_g303 r941_prom_fast_g507 r941_prom_hac_g507 r941_prom_high_g303 r941_prom_high_g330 r941_prom_high_g344 r941_prom_high_g360 r941_prom_high_g4011 r941_prom_sup_g507 r941_sup_plant_g610

  ```
  ### Command line example for creating a file inside a container:
  docker run --rm -v=/home/maksim/BIOINFOBSEE-461/Input_data/bacillus_assemblies_2022:/Input -v=/home/maksim/BIOINFOBSEE-461/Bt_19_BacPack:/Output -v=/home/maksim/BIOINFOBSEE-461/Ref:/Reference bacpack_renewed /Soft/BacPack/ConfigGen.sh --config /Output/Bt_19_bacpack_config.yaml --species_tag Bt_19 --fasta /Input/19.fna --output_dir /Output --ref_nucl_fasta_path /Reference --min_len 500
  ```
## Step 2: Running data analysis with a pipeline inside a container
 
  ```
  ### Command line example for running analysis inside a container:
  docker run --rm -v=/home/maksim/BIOINFOBSEE-461/Input_data/bacillus_assemblies_2022:/Input -v=/home/maksim/BIOINFOBSEE-461/Bt_19_BacPack:/Output -v=/home/maksim/BIOINFOBSEE-461/Ref:/Reference bacpack_renewed snakemake --snakefile /Soft/BacPack/BacPack.smk --cores 20 --configfile /Output/Bt_19_bacpack_config.yaml --use-conda
  ```
  
