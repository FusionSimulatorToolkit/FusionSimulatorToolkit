#!/bin/bash

# Template-based pass-dependent  PacBio read simulation using pbsim3
# Produces readset replicas for a range of: tissues, coverages, and number of passes
# Required software: pbsim3, samtools, pbindex, ccs, bedtools

simulateCoverage() {
    # creates a new FASTA file with multiple copies of each sequence
    local FASTA_IN=$1
    local COVERAGE=$2
    local FASTA_OUT=$3 
    echo --- Coverage simulation ---
    echo Input FATSA: ${FASTA_IN}
    echo Coverage: ${COVERAGE}
    echo Output FASTA: ${FASTA_OUT}
    if [[ -s ${FASTA_OUT}  ]]
    then 
        echo Output FATSA already exists: ${FASTA_OUT}
    else 
        for i in $(seq 1 $COVERAGE); do
            cat ${FASTA_IN} >> ${FASTA_OUT}
        done
    fi
}

simulatePacBioReads() {
  # simulates PacBio reads using PBSIM3* (*modified to report a name map)
  local FASTA=$1
  local TISSUE=$2
  local COVERAGE=$3
  local PASSES=$4
  local OUT_DIR=$5
  local PBSIM3_DIR=$6
  local PREFIX=${OUT_DIR}/${TISSUE}_cov${COVERAGE}_pass${PASSES}
  eval "$7='${PREFIX}.fastq.gz'"
  eval "$8='${PREFIX}.name_map'"
  echo --- PacBio Read Simulation: Start ---
  echo Input FATSA: ${FASTA} 
  echo Number of passes: ${PASSES}
  echo Output FASTQ: ${PREFIX}.fastq.gz  
  ${PBSIM3_DIR}/src/pbsim --strategy templ --method errhmm --errhmm ${PBSIM3_DIR}/data/ERRHMM-SEQUEL.model --template ${FASTA} --pass-num ${PASSES} --prefix ${PREFIX} --id-prefix ${TISSUE} &> ${PREFIX}.log
  
  if [ ${PASSES} -ne 1 ]; then
    samtools view -b ${PREFIX}.sam > ${PREFIX}.bam 
    rm ${PREFIX}.sam
    pbindex ${PREFIX}.bam
    ccs ${PREFIX}.bam ${PREFIX}.ccs.bam --min-passes ${PASSES} --min-rq 0
    bedtools bamtofastq -i ${PREFIX}.ccs.bam -fq ${PREFIX}.fastq
  fi
  gzip -f ${PREFIX}.fastq
  echo --- Read Simulation: Done ---
}


simulateONTReads() {
  # simulates ONT reads using the PBSIM3* R10 model, accuracy mean of 98% (*modified to report a name map)
  local FASTA=$1
  local TISSUE=$2
  local COVERAGE=$3
  local OUT_DIR=$4
  local PBSIM3_DIR=$5
  local PREFIX=${OUT_DIR}/${TISSUE}_cov${COVERAGE}
  eval "$6='${PREFIX}.fastq.gz'"
  eval "$7='${PREFIX}.name_map'"
  echo --- ONT Read Simulation: Start ---
  echo Input FATSA: ${FASTA} 
  echo Output FASTQ: ${PREFIX}.fastq.gz  
  ${PBSIM3_DIR}/src/pbsim --strategy templ --method errhmm --errhmm ${PBSIM3_DIR}/data/ERRHMM-ONT-HQ.model --template ${FASTA} --prefix ${PREFIX} --id-prefix ${TISSUE} --accuracy-mean 0.98 --pass 1 &> ${PREFIX}.log
  gzip -f ${PREFIX}.fastq
  echo --- Read Simulation: Done ---
}


#-------------------------#
#----------CONFIG---------#
#-------------------------#


usage () {
    echo USAGE: $0 input_fasta_dir output_dir pbsim_install_dir platform
    exit 1
}
if [ "$#" -ne 4 ]
then
    usage
fi

# directory with input FASTA files for each tissue
# expected name format: sim_${TISSUE}.fusion_transcripts.fasta
INPUT_DATA=$1 
# output directory where all the data will be written
OUTPUT_DIR=$2
# path to the PBSIM3 (updated) installation
PBSIM3_DIR=$3
# sequencing platform (PACBIO or ONT)
PLATFORM=$4
# number of replicas for each combination of tissue/coverage/passes
NREPLICAS=3
# tissue list (expected in the input directory)
TISSUES=('adipose' 'brain' 'colon' 'heart' 'testis')
# range of coverages to simulate
COVERAGES=(5 50)
# (PacBio only) number of passes to simulate at each coverage
if [ "$PLATFORM" = "PACBIO" ]; then
    PASSES=(1 2 3 4 5 10 20)
else
    PASSES=(1)
fi

#-------------------------#
#-------PIPELINE----------#
#-------------------------#

mkdir -p ${OUTPUT_DIR}/FASTA
mkdir -p ${OUTPUT_DIR}/READS
mkdir -p ${OUTPUT_DIR}/TISSUE_MIX

# clean up any files from a previous run
rm -rf ${OUTPUT_DIR}/READS/*
rm -rf ${OUTPUT_DIR}/TISSUE_MIX/*

echo -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# 1. prepare the FASTA files for each coverage
for COVERAGE in "${COVERAGES[@]}"; do
    for TISSUE in "${TISSUES[@]}"; do
        FASTA=${INPUT_DATA}/sim_${TISSUE}.fusion_transcripts.fasta
        FASTA_OUT=${OUTPUT_DIR}/FASTA/$(basename "${FASTA%.*}").cov${COVERAGE}.${FASTA##*.}
        simulateCoverage ${FASTA} ${COVERAGE} ${FASTA_OUT}
    done
done
echo -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# 2. generate reads
for COVERAGE in "${COVERAGES[@]}"; do
    for PASS in "${PASSES[@]}"; do
        for REPLICA in $(seq 1 $NREPLICAS); do
            mkdir -p ${OUTPUT_DIR}/READS/${REPLICA}
            mkdir -p ${OUTPUT_DIR}/TISSUE_MIX/${REPLICA}
	    for TISSUE in "${TISSUES[@]}"; do
                echo ----------------Coverage: ${COVERAGE} Number of passes: ${PASS} Replica: $REPLICA Tissue: $TISSUE
		FASTA=${OUTPUT_DIR}/FASTA/sim_${TISSUE}.fusion_transcripts.cov${COVERAGE}.fasta
		FASTQ=''
		NAME_MAP=''
                if [ "$PLATFORM" = "PACBIO" ]; then
		    simulatePacBioReads ${FASTA} ${TISSUE} ${COVERAGE} $PASS ${OUTPUT_DIR}/READS/${REPLICA} ${PBSIM3_DIR} FASTQ NAME_MAP 
	        else
                    simulateONTReads ${FASTA} ${TISSUE} ${COVERAGE} ${OUTPUT_DIR}/READS/${REPLICA} ${PBSIM3_DIR} FASTQ NAME_MAP
		fi
		# adds the reads to the tissue mix for this replica
		zcat ${FASTQ} >> ${OUTPUT_DIR}/TISSUE_MIX/${REPLICA}/mix_cov${COVERAGE}_pass${PASS}.fastq
                cat ${NAME_MAP} >> ${OUTPUT_DIR}/TISSUE_MIX/${REPLICA}/mix_cov${COVERAGE}_pass${PASS}_name_map.txt
            done
            gzip -f ${OUTPUT_DIR}/TISSUE_MIX/${REPLICA}/mix_cov${COVERAGE}_pass${PASS}.fastq
        done
    done
done
echo -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
