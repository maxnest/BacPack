#!/bin/bash

echo "                                                                                                                                   
                                                                                                                                   
BBBBBBBBBBBBBBBBB                                        PPPPPPPPPPPPPPPPP                                      kkkkkkkk           
B::::::::::::::::B                                       P::::::::::::::::P                                     k::::::k           
B::::::BBBBBB:::::B                                      P::::::PPPPPP:::::P                                    k::::::k           
BB:::::B     B:::::B                                     PP:::::P     P:::::P                                   k::::::k           
  B::::B     B:::::B  aaaaaaaaaaaaa      cccccccccccccccc  P::::P     P:::::Paaaaaaaaaaaaa      cccccccccccccccc k:::::k    kkkkkkk
  B::::B     B:::::B  a::::::::::::a   cc:::::::::::::::c  P::::P     P:::::Pa::::::::::::a   cc:::::::::::::::c k:::::k   k:::::k 
  B::::BBBBBB:::::B   aaaaaaaaa:::::a c:::::::::::::::::c  P::::PPPPPP:::::P aaaaaaaaa:::::a c:::::::::::::::::c k:::::k  k:::::k  
  B:::::::::::::BB             a::::ac:::::::cccccc:::::c  P:::::::::::::PP           a::::ac:::::::cccccc:::::c k:::::k k:::::k   
  B::::BBBBBB:::::B     aaaaaaa:::::ac::::::c     ccccccc  P::::PPPPPPPPP      aaaaaaa:::::ac::::::c     ccccccc k::::::k:::::k    
  B::::B     B:::::B  aa::::::::::::ac:::::c               P::::P            aa::::::::::::ac:::::c              k:::::::::::k     
  B::::B     B:::::B a::::aaaa::::::ac:::::c               P::::P           a::::aaaa::::::ac:::::c              k:::::::::::k     
  B::::B     B:::::Ba::::a    a:::::ac::::::c     ccccccc  P::::P          a::::a    a:::::ac::::::c     ccccccc k::::::k:::::k    
BB:::::BBBBBB::::::Ba::::a    a:::::ac:::::::cccccc:::::cPP::::::PP        a::::a    a:::::ac:::::::cccccc:::::ck::::::k k:::::k   
B:::::::::::::::::B a:::::aaaa::::::a c:::::::::::::::::cP::::::::P        a:::::aaaa::::::a c:::::::::::::::::ck::::::k  k:::::k  
B::::::::::::::::B   a::::::::::aa:::a cc:::::::::::::::cP::::::::P         a::::::::::aa:::a cc:::::::::::::::ck::::::k   k:::::k 
BBBBBBBBBBBBBBBBB     aaaaaaaaaa  aaaa   ccccccccccccccccPPPPPPPPPP          aaaaaaaaaa  aaaa   cccccccccccccccckkkkkkkk    kkkkkkk
                                                                                                                                   
                                                                                                                                  "

usage() { echo "Usage: $0 --config <string> --species_tag <string> --r1_fq <string> --r2_fq <string> --output_dir <string>
       --ref_nucl_fasta_path <string> --min_len <string> [--long_fq <string>] [--medaka_model <string>] " 1>&2; }

eval set -- `getopt --options '' --longoptions config:,species_tag:,r1_fq:,r2_fq:,output_dir:,ref_nucl_fasta_path:,min_len:,long_fq::,medaka_model:: -- $@`

while true; do
    case ${1} in
	--config)
	    config=${2}
	    echo "*** All arguments will be written to a file ${config} ***"
	    shift 2 ;;
        --species_tag)
	    species_tag=${2}
	    echo "*** The ${species_tag} is set as the value of the 'species_tag' argument ***"
	    shift 2 ;;
        --r1_fq)
	    r1_fq=${2}
	    echo "*** The ${r1_fq} is set as the short paired-end forward (R1) library ***"
	    shift 2 ;;
        --r2_fq)
	    r2_fq=${2}
	    echo "*** The ${r2_fq} is set as the short paired-end reverse (R2) library ***"
	    shift 2 ;;
        --output_dir)
	    output_dir=${2}
	    echo "*** Directory ${output_dir} selected as output directory for results obtained ***"
	    shift 2 ;;
	--ref_nucl_fasta_path)
	    ref_nucl_fasta_path=${2}
	    echo "*** Directory ${ref_nucl_fasta_path} selected as containing reference nucleotide sequences ***"
	    shift 2 ;;
	--min_len)
            min_len=${2}
	    echo "*** ${min_len} nucleotides set as the minimum length of the assembled sequence ***"
	    shift 2 ;;
	--long_fq)
	    long_fq=${2}
	    shift 2 ;;
	--medaka_model)
	    medaka_model=${2}
	    shift 2 ;;
	--)
	    break
	    shift ;;
	*)
	    usage
	    exit 1
	    ;;
    esac
done

## MESSAGE 
if [ -z ${config} ] && [ -z ${species_tag} ] && [ -z ${r1_fq} ] && [ -z ${r2_fq} ] && [ -z ${output_dir} ] && [ -z ${ref_nucl_fasta_path} ] && [ -z ${min_len} ] && [ -z ${long_fq} ] && [ -z ${medaka_model} ] ; then
	usage
        echo "
	--config	The name of the configuration file in which the parameters for launching the pipeline will be written. As an example, 'BacPack_config.yaml'
	--species_tag	Species and/or strain name to be used when generating analysis results. As an example, 'Bt_st1402'
	--r1_fq	Full path to the archive with raw forward (R1) short paired-end reads
	--r2_fq Full path to the archive with raw reverse (R2) short paired-end reads 
	--output_dir The name of the directory where all the results of the pipeline will be written
	--ref_nucl_fasta_path	Full path to the directory containing fasta-files with genome assemblies that will be considered as references
	--min_len	Minimum length of assembled sequences
	--long_fq	Full path to the archive with raw long reads
	--medaka_model	The name of the model that will be used to create consensus based on the results of the assembly of long reads. Available options: r103_fast_g507 r103_hac_g507 r103_min_high_g345 r103_min_high_g360 r103_prom_high_g360 r103_sup_g507 r1041_e82_260bps_fast_g632 r1041_e82_260bps_hac_g632 r1041_e82_260bps_sup_g632 r1041_e82_400bps_fast_g615 r1041_e82_400bps_fast_g632 r1041_e82_400bps_hac_g615 r1041_e82_400bps_hac_g632 r1041_e82_400bps_sup_g615 r104_e81_fast_g5015 r104_e81_hac_g5015 r104_e81_sup_g5015 r104_e81_sup_g610 r10_min_high_g303 r10_min_high_g340 r941_e81_fast_g514 r941_e81_hac_g514 r941_e81_sup_g514 r941_min_fast_g303 r941_min_fast_g507 r941_min_hac_g507 r941_min_high_g303 r941_min_high_g330 r941_min_high_g340_rle r941_min_high_g344 r941_min_high_g351 r941_min_high_g360 r941_min_sup_g507 r941_prom_fast_g303 r941_prom_fast_g507 r941_prom_hac_g507 r941_prom_high_g303 r941_prom_high_g330 r941_prom_high_g344 r941_prom_high_g360 r941_prom_high_g4011 r941_prom_sup_g507 r941_sup_plant_g610"
	exit 1
fi

## WRITING ARGUMENTS TO A CONFIGURATION FILE
echo "# Input: #" > ${config}
 
if [ ! -z ${species_tag} ]; then
    echo "species_tag: ${species_tag}" >> ${config}
fi

if [ ! -z ${r1_fq} ] && [ ! -z ${r2_fq} ]; then
    echo "r1_fq: ${r1_fq}" >> ${config}
    echo "r2_fq: ${r2_fq}" >> ${config}
fi

if [ ! -z ${long_fq} ]; then
    echo "*** The ${long_fq} is set as the long read library ***"
    echo "long_fq: ${long_fq}" >> ${config}
else
    echo "long_fq: FALSE" >> ${config}
fi

if [ ! -z ${medaka_model} ]; then
    echo "*** The ${medaka_model} model is selected as the model to run the Medaka ***"
    echo "medaka_model: ${medaka_model}" >> ${config}
else
    echo "medaka_model: FALSE" >> ${config}
fi

if [ ! -z ${output_dir} ]; then
    echo "output_dir: ${output_dir}" >> ${config}
fi

if [ ! -z ${ref_nucl_fasta_path} ]; then
    echo "ref_nucl_fasta_path: ${ref_nucl_fasta_path}" >> ${config}
fi

if [ ! -z ${min_len} ]; then
    echo "min_len: ${min_len}" >> ${config}
fi

# RESOURCES
echo "# Resources: #" >> ${config}
echo "threads: 20" >> ${config}
echo "prokka_db: /Soft/BacPack/resources/IPG_Bt_sequence.fasta" >> ${config}
#echo "bacillales_odb: /home/maksim/Soft/BUSCO_db/bacillales_odb10" >> ${config}
#echo "bacilli_odb: /home/maksim/Soft/BUSCO_db/bacilli_odb10" >> ${config}

# Soft
echo "# Soft: #" >> ${config}
echo "python: /opt/conda/bin/python" >> ${config}
echo "fastqc: /Soft/FastQC/fastqc" >> ${config}
echo "fastp:  /Soft/fastp" >> ${config}
echo "spades: /Soft/SPAdes-3.15.4-Linux/bin/spades.py" >> ${config}
echo "quast: /Soft/quast-5.2.0/quast.py" >> ${config}
#echo "busco: /home/maksim/anaconda3/envs/busco_env/bin/busco" >> ${config}
echo "fastani: /Soft/fastANI" >> ${config}
echo "prodigal: /Soft/prodigal.linux" >> ${config}
echo "prokka: /Soft/prokka/bin/prokka" >> ${config}
echo "checkm: /opt/conda/envs/checkm/bin/checkm" >> ${config}
echo "cryprocessor: /Soft/cry_processor/cry_processor.py" >> ${config}
#echo "idops: /home/maksim/anaconda3/envs/idops/bin/idops" >> ${config}
echo "rabbitqc: /Soft/RabbitQCPlus/RabbitQCPlus" >> ${config}
echo "flye: /Soft/Flye/bin/flye" >> ${config}
echo "medaka_consensus: /Soft/medaka/venv/bin/medaka_consensus" >> ${config}
echo "bwa: /Soft/bwa/bwa" >> ${config}
echo "samtools: /Soft/samtools-1.16.1/samtools" >> ${config}
echo "pilon: /Soft/pilon-1.24.jar" >> ${config}
echo "deepbgc: /opt/conda/envs/deepbgc/bin/deepbgc" >> ${config}
echo "/opt/conda/envs/antismash_env/bin/antismash" >> ${config}
