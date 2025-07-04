#!/bin/bash

# Define input and output files
input_file="afr_mtag.tsv"
output_file="fuma_jass_afr.txt"
variantqc_file="afr_eur_variant_qc_metrics.txt "

# Write header to the output file
echo "SNP,Zscore,P_value,freq,SE,beta,varbeta" > $output_file

# Read the input file line by line (skipping the header)
tail -n +2 $input_file | while IFS=$'\t' read -r Region CHR snp_ids position Ref_allele Alt_allele MiddlePosition JASS_PVAL UNIVARIATE_MIN_PVAL UNIVARIATE_MIN_QVAL PLEIOTROPY_INDEX z_UKB_PULSERATE z_UKB_MAXHEARTRATE z_UKB_SYSTOLICBP z_UKB_DIASTOLICBP
do
    # Create SNP field
    SNP="${CHR}:${position}_${Ref_allele}_${Alt_allele}"
    
    # Zscore is calculated by averaging the z-scores from different traits
    Zscore=$(echo "($z_UKB_PULSERATE + $z_UKB_MAXHEARTRATE + $z_UKB_SYSTOLICBP + $z_UKB_DIASTOLICBP) / 4" | bc -l)

    # Use JASS_PVAL for P_value
    P_value=$JASS_PVAL

    # Get frequency from variant QC file
    freq=$(awk -v snp_id=$snp_ids 'BEGIN{FS="\t";OFS="\t"} $1==snp_id {print $3}' $variantqc_file)

    # Calculate SE using the formula from MATLAB code
    SE=$(echo "scale=6; $Zscore / (p_value/2)" | bc -l)

    # Estimate beta using Zscore and SE
    beta=$(echo "scale=6; $Zscore * $SE" | bc -l)

    # Estimate varbeta using SE
    varbeta=$(echo "scale=6; $SE * $SE" | bc -l)

    # Write the formatted line to the output file
    echo "$SNP,$Zscore,$P_value,$freq,$SE,$beta,$varbeta" >> $output_file
done
