#!/bin/bash
# title: eQTL pipeline
# authors: Fabio Marroni, Andrea Bertoli, Giovanni Gabelli
# date: 04/06/2018


############################################################################################################
#PRELIMINARY STEPS
############################################################################################################
READSPATH=/iga/illumina-spool/FC0903/reads/uniud-novabreed                                      #folder with original reads
EQTLSEQDIR=/projects/novabreed/lanes/sequences/rna/vitis_vinifera/eQTL_var                      #folder with copied files and trimmed files
EQTLALIGN=/projects/novabreed/lanes/alignments/STAR/rna/vitis_vinifera/genome_12xCHR/eQTL_var   #folder with STAR output
EQTLRESULTS=/projects/novabreed/share/abertoli/eQTL/files2                                      #folder with the results
RPATH=/projects/novabreed/share/abertoli/eQTL/script/R                                          #folder with the R scripts
#put SNP files in $EQTLRESULTS, put SNP and genes position files in $EQTLRESULTS

mkdir $EQTLRESULTS/1_upstream_matrixEQTL
mkdir $EQTLRESULTS/2_input_matrixEQTL
mkdir $EQTLRESULTS/3_output_matrixEQTL
mkdir $EQTLRESULTS/4_downstream_matrixEQTL
mkdir $EQTLRESULTS/5_GOterm_extraction
mkdir $EQTLRESULTS/function_annot               #put the GO-gene file here
mkdir $EQTLRESULTS/herit
mkdir $EQTLRESULTS/star_counts                  #put STAR gene counts here
mkdir $EQTLRESULTS/glm                          #put glm table here

############################  0. Create a report with number of reads after trimming, after filtering and after alignment ###########################
REPORT_PATH=/projects/novabreed/share/abertoli/eQTL/files2/1_upstream_matrixEQTL
rm $REPORT_PATH/trimFilt_report.txt;

echo "Trimming report, only paired reads are considered. Numbers indicate a single read (not a pair) that HAVE a mate." $'\n' "SAMPLE_id" $'\t' "#READS_after_trim" $'\t' "#READS_after_filt" $'\t' "#READS_after_alignment" >$REPORT_PATH/trimFilt_report.txt

for FULLPATH in $EQTLSEQDIR/*
do
  SAMPLE=$(basename $FULLPATH)
  echo "Working on: $SAMPLE"

  for R in 1 2   #per ogni replica
    do
    REP="rep$R"
    nREP=$(ls ${FULLPATH}/*$REP*R1.fastq.gz | head -n1)
    if [ -f $nREP  ]                                          #if $nREP exists
      then
        R_BEFORE_FILT=$(( $(zcat $FULLPATH/*$REP*R1* | wc -l) / 2 ))
        R_AFTER_FILT=$(( $(zcat $FULLPATH/trimmed/*$REP*_1*.fastq.gz | wc -l) / 2 ))

        nADD=$(ls $EQTLALIGN/$SAMPLE/*${REP}*final.out | wc -l)
        if [ $nADD -eq 1 ]
          then
            R_UNIQUE_ALIGN=$(( $(grep "Uniquely mapped reads number" $EQTLALIGN/$SAMPLE/*$REP*final.out  | cut -d"|" -f2) * 2 ))
          fi
        if [ $nADD -eq 2 ]
          then
            ADD1=$(echo $(grep "Uniquely mapped reads number" $EQTLALIGN/$SAMPLE/*$REP*final.out  | cut -d"|" -f2 | head -n1))
            ADD2=$(echo $(grep "Uniquely mapped reads number" $EQTLALIGN/$SAMPLE/*$REP*final.out  | cut -d"|" -f2 | tail -n1))
            R_UNIQUE_ALIGN=$(( $(($ADD1 + $ADD2)) * 2))
          fi

        touch $REPORT_PATH/tmp1.txt
        touch $REPORT_PATH/tmp2.txt
        echo "${SAMPLE}_${REP}" $'\t' $R_BEFORE_FILT $'\t' $R_AFTER_FILT $'\t' $R_UNIQUE_ALIGN > $REPORT_PATH/tmp1.txt
        cat $REPORT_PATH/trimFilt_report.txt $REPORT_PATH/tmp1.txt > $REPORT_PATH/tmp2.txt
        mv $REPORT_PATH/tmp2.txt $REPORT_PATH/trimFilt_report.txt
      fi
  done

done
rm $REPORT_PATH/*tmp*

############################  1. Merge multiples run from the same sample, put those files in a subdirectory ###########################
mkdir $EQTLRESULTS/star_counts/merged/
for R in 1 2   #for each replicate
  do
  REP="rep$R"
  for DIRECTORY in $EQTLALIGN/*  #for each cultivar
    do
    NAME=$(basename $DIRECTORY)
    NEW_NAME=$(echo "$NAME" | tr _ -)
    echo $NEW_NAME
    x=$(ls $EQTLRESULTS/star_counts/${NEW_NAME}_${REP}* | wc -l)
    if [ $x -eq 1 ]       #se non ci sono aggiunte per questa replica
      then
      y=$(ls $EQTLRESULTS/star_counts/${NEW_NAME}_${REP}*)
      echo $y
      cp $y $EQTLRESULTS/star_counts/merged/${NEW_NAME}_${REP}.txt
      fi
    if [ $x -eq 2 ]       #se ci sono aggiunte per questa replica
      then
      touch $EQTLRESULTS/star_counts/tmp1.txt
      touch $EQTLRESULTS/star_counts/tmp2.txt
      touch $EQTLRESULTS/star_counts/tmp3.txt
      for TAB in $EQTLRESULTS/star_counts/${NEW_NAME}_${REP}*     #per ogni file della replica
      do
        echo $TAB
        cut -f 2 $TAB > $EQTLRESULTS/star_counts/tmp1.txt
        paste $EQTLRESULTS/star_counts/tmp1.txt $EQTLRESULTS/star_counts/tmp2.txt > $EQTLRESULTS/star_counts/tmp3.txt
        cut -f 1,2 $EQTLRESULTS/star_counts/tmp3.txt > $EQTLRESULTS/star_counts/tmp2.txt
        rm $EQTLRESULTS/star_counts/tmp3.txt
      done
      sed 's/\t/+/g' $EQTLRESULTS/star_counts/tmp2.txt | bc > $EQTLRESULTS/star_counts/tmp3.txt
      cut -f 1 $(ls $EQTLRESULTS/star_counts/${NEW_NAME}_${REP}* | head -1) >  $EQTLRESULTS/star_counts/gnames.txt
      paste  $EQTLRESULTS/star_counts/gnames.txt $EQTLRESULTS/star_counts/tmp3.txt > $EQTLRESULTS/star_counts/merged/${NEW_NAME}_${REP}.txt
      rm $EQTLRESULTS/star_counts/gnames.txt
      rm $EQTLRESULTS/star_counts/tmp*
      fi
    if [ $x -eq 0 ] #la replica per questa varietà non esiste
      then
        echo "${REP} not present for ${NEW_NAME}"
      fi
    echo ""
  done
done


############################################################################################################
#RNA-SEQ DATA ANALYSIS
############################################################################################################
#detected wrong sample, we have to delete them
rm $EQTLRESULTS/star_counts/merged/lambrusco-grasparossa_rep2*
rm $EQTLRESULTS/star_counts/merged/marandi*

###########################  2. Expression data normalization  ############################
IN_DIR=$EQTLRESULTS/star_counts/merged
OUT_DATA=$EQTLRESULTS/expN_data.txt
R --no-save --args a1=$IN_DIR a2=$OUT_DATA < $RPATH/A_exp_normalization.R

###########################  3. Varieties clustering (distance=corr.Spearman, method=complete-linkage) ##########################
IN_DATA=$EQTLRESULTS/expN_data.txt
OUT_DATA=$EQTLRESULTS/1_upstream_matrixEQTL/clustering.pdf
R --no-save --args a1=$IN_DATA a2=$OUT_DATA < $RPATH/B_clustering.R

###########################  4. PCA with genotypes and expression data  ############################
INPUT_SNP=$EQTLRESULTS/snp_data.txt
INPUT_EXP=$EQTLRESULTS/expN_data.txt
OUTPUT_FOLDER=$EQTLRESULTS/1_upstream_matrixEQTL/
R --no-save --args a1=$INPUT_SNP a2=$INPUT_EXP a3=$OUTPUT_FOLDER < $RPATH/C_PCA.R

###########################  5. Remove wrong replicate and merge good replicates  ###########################
REMOVED_REP=montepulciano_rep1
INPUT_FILE=$EQTLRESULTS/expN_data.txt
OUTPUT_FILE=$EQTLRESULTS/expN_data_merged.txt
R --no-save --args a1=$REMOVED_REP a2=$INPUT_FILE a3=$OUTPUT_FILE < $RPATH/D_exp_merging.R

###########################  6. Expression data filtering (+summary)  ###########################
INPUT_FILE=$EQTLRESULTS/expN_data_merged.txt
OUTPUT_FILE=$EQTLRESULTS/expN_data_merged_filtered.txt
VAR_MIN=0.1
MEDIAN_MIN=1
R --no-save --args a1=$INPUT_FILE a2=$OUTPUT_FILE a3=$VAR_MIN a4=$MEDIAN_MIN < $RPATH/E_exp_filtering.R

###########################  7. Heritability  ###########################
INPUT_FILE=$EQTLRESULTS/expN_data.txt
FILTERED_DATA=$EQTLRESULTS/expN_data_merged_filtered.txt
REMOVED_REP=montepulciano_rep1
OUTPUT_FOLDER=$EQTLRESULTS/herit/
R --no-save --args a1=$INPUT_FILE a2=$FILTERED_DATA a3=$REMOVED_REP a4=$OUTPUT_FOLDER < $RPATH/F_heritability.R


############################################################################################################
#SNPs DATA ANALYSIS
############################################################################################################

########################### 8. SNPs data filtering (+summary) ###########################
INPUT_FILE=$EQTLRESULTS/snp_data.txt
OUTPUT_FILE=$EQTLRESULTS/snp_data_filtered.txt
MIN_nGT=0.3
MIN_GF=8
REMOVED_VARIETIES=marandi.shemakhinskii.564         #for SNPs
R --no-save --args a1=$INPUT_FILE a2=$OUTPUT_FILE a3=$MIN_nGT a4=$MIN_GF a5=$REMOVED_VARIETIES < $RPATH/G_SNP_filtering.R


############################################################################################################
#MATRIX eQTL ANALYSIS
############################################################################################################

########################### 9. Prepare input files for matrix eQTL analysis ###########################
SNP_DATA=$EQTLRESULTS/snp_data_filtered.txt
EXP_DATA=$EQTLRESULTS/expN_data_merged_filtered.txt
FORMATTED_SNP=$EQTLRESULTS/2_input_matrixEQTL/snp_data_formatted.txt
FORMATTED_EXP=$EQTLRESULTS/2_input_matrixEQTL/exp_data_formatted.txt
GENE_POS=$EQTLRESULTS/gene_pos.txt
SNP_POS=$EQTLRESULTS/snp_pos.txt
FORMATTED_GENE_POS=$EQTLRESULTS/2_input_matrixEQTL/gene_pos.txt
FORMATTED_SNP_POS=$EQTLRESULTS/2_input_matrixEQTL/snp_pos.txt
COVAR=$EQTLRESULTS/1_upstream_matrixEQTL/PC_snp.txt
FORMATTED_COVAR=$EQTLRESULTS/2_input_matrixEQTL/PC_SNP_mEQTL.txt
N_PC=5
REMOVED_VARITIES=marandi.shemakhinskii.564         #for CVs
R --no-save --args a1=$SNP_DATA a2=$EXP_DATA a3=$FORMATTED_SNP a4=$FORMATTED_EXP a5=$GENE_POS a6=$SNP_POS a7=$FORMATTED_GENE_POS a8=$FORMATTED_SNP_POS a9=$COVAR a10=$FORMATTED_COVAR a11=$N_PC a12=$REMOVED_VARITIES < $RPATH/H_matrixEQTL_prep.R

########################### 10. Matrix eQTL analysis ###########################
FORMATTED_SNP=$EQTLRESULTS/2_input_matrixEQTL/snp_data_formatted.txt
FORMATTED_EXP=$EQTLRESULTS/2_input_matrixEQTL/exp_data_formatted.txt
FORMATTED_GENE_POS=$EQTLRESULTS/2_input_matrixEQTL/gene_pos.txt
FORMATTED_SNP_POS=$EQTLRESULTS/2_input_matrixEQTL/snp_pos.txt
FORMATTED_COVAR=$EQTLRESULTS/2_input_matrixEQTL/PC_SNP_mEQTL.txt
OUTPUT_FOLDER=$EQTLRESULTS/3_output_matrixEQTL/
P_VALUE=1e-7

R --no-save --args a1=$SNP_DATA a2=$EXP_DATA a3=$FORMATTED_GENE_POS a4=$FORMATTED_SNP_POS a5=$FORMATTED_COVAR a6=$OUTPUT_FOLDER a7=$P_VALUE < $RPATH/I_matrixEQTL_analysis.R

########################### 11a. Boxplot for eGene ###########################
GENE_TO_CHECK=VIT_206s0004g06490
SNP_FILE=$EQTLRESULTS/snp_data_filtered.txt
EXP_FILE=$EQTLRESULTS/expN_data_merged_filtered.txt
RESULTS_FOLDER=$EQTLRESULTS/3_output_matrixEQTL/
OUTPUT_FOLDER=$EQTLRESULTS/4_downstream_matrixEQTL/
R --no-save --args a1=$GENE_TO_CHECK a2=$SNP_FILE a3=$EXP_FILE a4=$RESULTS_FOLDER a5=$OUTPUT_FOLDER < $RPATH/Z1_boxplot_gene.R

########################### 11b. Boxplot for eQTL ###########################
SNP_TO_CHECK=chr4_21392549
SNP_FILE=$EQTLRESULTS/snp_data_filtered.txt
EXP_FILE=$EQTLRESULTS/expN_data_merged_filtered.txt
RESULTS_FOLDER=$EQTLRESULTS/3_output_matrixEQTL/
OUTPUT_FOLDER=$EQTLRESULTS/4_downstream_matrixEQTL/
R --no-save --args a1=$SNP_TO_CHECK a2=$SNP_FILE a3=$EXP_FILE a4=$RESULTS_FOLDER a5=$OUTPUT_FOLDER < $RPATH/Z2_boxplot_SNP.R

########################### 12. eQTL results analysis: final plot, GO terms ###########################
#final plot for all eQTL results
RESULTS_TRANS=$EQTLRESULTS/3_output_matrixEQTL/results_trans.txt
RESULTS_CIS=$EQTLRESULTS/3_output_matrixEQTL/results_cis.txt
GENE_POS=$EQTLRESULTS/2_input_matrixEQTL/gene_pos.txt
SNP_POS=$EQTLRESULTS/2_input_matrixEQTL/snp_pos.txt
OUTPUT_FILE=$EQTLRESULTS/3_output_matrixEQTL/eQTL_results_plot.png
ACTION=normal
R --no-save --args a1=$RESULTS_TRANS a2=$RESULTS_CIS a3=$GENE_POS a4=$SNP_POS a5=$OUTPUT_FILE a6=NULL a7=$ACTION < $RPATH/L_final_plot.R

#GO terms analysis for cis results
GENE_TO_GO=$EQTLRESULTS/function_annot/GO.txt
INPUT_INT_GENES=$EQTLRESULTS/3_output_matrixEQTL/results_cis.txt
INPUT_ALL_GENES=$EQTLRESULTS/expN_data_merged_filtered.txt
OUTPUT_PREFIX=$EQTLRESULTS/3_output_matrixEQTL/GO_cis
FIRST_GO=15
R --no-save --args a1=$GENE_TO_GO a2=$INPUT_INT_GENES a3=$INPUT_ALL_GENES a4=$OUTPUT_PREFIX a5=$FIRST_GO < $RPATH/M_GOtest_eQTLresults.R

#GO terms analysis for trans results
GENE_TO_GO=$EQTLRESULTS/function_annot/GO.txt
INPUT_INT_GENES=$EQTLRESULTS/3_output_matrixEQTL/results_trans.txt
INPUT_ALL_GENES=$EQTLRESULTS/expN_data_merged_filtered.txt
OUTPUT_PREFIX=$EQTLRESULTS/3_output_matrixEQTL/GO_trans
FIRST_GO=15
R --no-save --args a1=$GENE_TO_GO a2=$INPUT_INT_GENES a3=$INPUT_ALL_GENES a4=$OUTPUT_PREFIX a5=$FIRST_GO <$RPATH/M_GOtest_eQTLresults.R


############################################################################################################
#IN DETAIL ANALYSIS OF PATHWAYS AND GENOMIC REGIONS
############################################################################################################

########################### 13. Specific pathway analysis ###########################
#GO term genes extraction and single plots
GO_TERM=GO:0009813
GO_TERM_OUTPUT=GO_0009813
GO_ANNOTATION=$EQTLRESULTS/function_annot/GO.txt
RESULT_CIS=$EQTLRESULTS/3_output_matrixEQTL/results_cis.txt
RESULT_TRANS=$EQTLRESULTS/3_output_matrixEQTL/results_trans.txt
OUTPUT_FOLDER=$EQTLRESULTS/5_GOterm_extraction/$GO_TERM_OUTPUT
GENE_POS=$EQTLRESULTS/2_input_matrixEQTL/gene_pos.txt
SNP_POS=$EQTLRESULTS/2_input_matrixEQTL/snp_pos.txt
mkdir $OUTPUT_FOLDER
R --no-save --args a1=$GO_TERM a2=$GO_ANNOTATION a3=$RESULT_CIS a4=$RESULT_TRANS a5=$OUTPUT_FOLDER a6=$GENE_POS a7=$SNP_POS a8=$GO_TERM_OUTPUT < $RPATH/N_GO_extraction.R

#big plot for the selected GO term
RESULTS_TRANS=$OUTPUT_FOLDER/${GO_TERM_OUTPUT}_trans.txt
RESULTS_CIS=$OUTPUT_FOLDER/${GO_TERM_OUTPUT}_cis.txt
GENE_POS=$EQTLRESULTS/2_input_matrixEQTL/gene_pos.txt
SNP_POS=$EQTLRESULTS/2_input_matrixEQTL/snp_pos.txt
OUTPUT_FILE=$OUTPUT_FOLDER/${GO_TERM_OUTPUT}_plot.png
ACTION=normal
R --no-save --args a1=$RESULTS_TRANS a2=$RESULTS_CIS a3=$GENE_POS a4=$SNP_POS a5=$OUTPUT_FILE a6=NULL a7=$ACTION < $RPATH/L_final_plot.R

########################### 14. Plot for a specific gene in a specific genome window ###########################
GENE=VIT_209s0070g00240
CHR=9
CENTER=15000000
RANGE=300000

RESULT_CIS=$EQTLRESULTS/3_output_matrixEQTL/results_cis.txt
RESULT_TRANS=$EQTLRESULTS/3_output_matrixEQTL/results_trans.txt
OUTPUT_FOLDER=$EQTLRESULTS/4_downstream_matrixEQTL
GENE_POS=$EQTLRESULTS/2_input_matrixEQTL/gene_pos.txt
SNP_POS=$EQTLRESULTS/2_input_matrixEQTL/snp_pos.txt
R --no-save --args a1=$GENE a2=$CHR a3=$CENTER a4=$RANGE a5=$RESULT_CIS a6=$RESULT_TRANS a7=$OUTPUT_FOLDER a8=$GENE_POS a9=$SNP_POS < $RPATH/O_p-value_plot.R


############################################################################################################
#CO-EXPRESSION NETWORK
############################################################################################################
mkdir $EQTLRESULTS/4_downstream_matrixEQTL/co-expression_network
mkdir $EQTLRESULTS/4_downstream_matrixEQTL/co-expression_network_analysis

########################### 15. Co-expression network construction ###########################
EXP_DATA=$EQTLRESULTS/expN_data_merged_filtered.txt
OUTPUT_FOLDER=$EQTLRESULTS/4_downstream_matrixEQTL/co-expression_network/
R --no-save --args a1=$EXP_DATA a2=$OUTPUT_FOLDER < $RPATH/P_Network.R

########################### 16. Co-expression network analysis (GO terms+summary) ###########################
for FILE in $EQTLRESULTS/4_downstream_matrixEQTL/co-expression_network/WGCNAnodes_*
do
  MODULE=$(basename $FILE)
  echo Starting analysis: $MODULE
  LIST_PATH=$EQTLRESULTS/4_downstream_matrixEQTL/co-expression_network_analysis/${MODULE//.txt/_list.txt}
  cat $FILE | tail -n +2 | cut -f1 > $LIST_PATH

  GENE_TO_GO=$EQTLRESULTS/function_annot/GO.txt
  INPUT_INT_GENES=$LIST_PATH
  INPUT_ALL_GENES=$EQTLRESULTS/expN_data_merged_filtered.txt
  OUTPUT_PREFIX=$EQTLRESULTS/4_downstream_matrixEQTL/co-expression_network_analysis/${MODULE//.txt/}
  FIRST_GO=15
  R --no-save --args a1=$GENE_TO_GO a2=$INPUT_INT_GENES a3=$INPUT_ALL_GENES a4=$OUTPUT_PREFIX a5=$FIRST_GO < $RPATH/Z3_GOtest.R
done

#generate summaries for each module + general summary
ALL_SUMMARY_FILE=$EQTLRESULTS/4_downstream_matrixEQTL/co-expression_network_analysis/modules_summaryALL.txt
rm $ALL_SUMMARY_FILE
TABLE_FILE=$EQTLRESULTS/4_downstream_matrixEQTL/co-expression_network_analysis/modules_tableALL.txt
rm $TABLE_FILE

for FILE in $EQTLRESULTS/4_downstream_matrixEQTL/co-expression_network/WGCNAnodes_*
do
  MODULE=$(basename $FILE)
  echo Summary \for: $MODULE
  LIST_PATH=$EQTLRESULTS/4_downstream_matrixEQTL/co-expression_network_analysis/${MODULE//.txt/_list.txt}
  N_GENES=$(cat $LIST_PATH | wc -l)
  echo ${MODULE//.txt/} $'\t' $N_GENES >>$TABLE_FILE

  LIST_GO_BP=$EQTLRESULTS/4_downstream_matrixEQTL/co-expression_network_analysis/${MODULE//.txt/_GO_resultsBP.txt}
  LIST_GO_MF=$EQTLRESULTS/4_downstream_matrixEQTL/co-expression_network_analysis/${MODULE//.txt/_GO_resultsMF.txt}
  LIST_GO_CC=$EQTLRESULTS/4_downstream_matrixEQTL/co-expression_network_analysis/${MODULE//.txt/_GO_resultsCC.txt}
  SUMMARY_FILE=$EQTLRESULTS/4_downstream_matrixEQTL/co-expression_network_analysis/${MODULE//.txt/_summary.txt}
  rm $SUMMARY_FILE

  echo "Number of genes": $N_GENES $'\n' >> $SUMMARY_FILE
  echo Biological process >> $SUMMARY_FILE
  cat $LIST_GO_BP >> $SUMMARY_FILE
  echo >> $SUMMARY_FILE
  echo Molecular Function >> $SUMMARY_FILE
  cat $LIST_GO_MF >> $SUMMARY_FILE
  echo >> $SUMMARY_FILE
  echo Cellular Component >> $SUMMARY_FILE
  cat $LIST_GO_CC >> $SUMMARY_FILE
  echo >> $SUMMARY_FILE

  echo ${MODULE//.txt/} - "Number of genes": $N_GENES $'\n' >> $ALL_SUMMARY_FILE
  echo Biological process >> $ALL_SUMMARY_FILE
  cat $LIST_GO_BP | head -n5 >> $ALL_SUMMARY_FILE
  echo >> $ALL_SUMMARY_FILE
  echo Molecular Function >> $ALL_SUMMARY_FILE
  cat $LIST_GO_MF | head -n5 >> $ALL_SUMMARY_FILE
  echo >> $ALL_SUMMARY_FILE
  echo Cellular Component >> $ALL_SUMMARY_FILE
  cat $LIST_GO_CC | head -n5 >> $ALL_SUMMARY_FILE
  echo >> $ALL_SUMMARY_FILE

done

########################### 17. Big plot + module colors ###########################
#print plot with all eQTL results, colors=modules
RESULTS_TRANS=$EQTLRESULTS/3_output_matrixEQTL/results_trans.txt
RESULTS_CIS=$EQTLRESULTS/3_output_matrixEQTL/results_cis.txt
GENE_POS=$EQTLRESULTS/2_input_matrixEQTL/gene_pos.txt
SNP_POS=$EQTLRESULTS/2_input_matrixEQTL/snp_pos.txt
OUTPUT_FILE=$EQTLRESULTS/4_downstream_matrixEQTL/eQTL+modules_results_plot.png
NET_FILE=$EQTLRESULTS/4_downstream_matrixEQTL/co-expression_network/WGCNAnodes.txt
ACTION=moduleColors
R --no-save --args a1=$RESULTS_TRANS a2=$RESULTS_CIS a3=$GENE_POS a4=$SNP_POS a5=$OUTPUT_FILE a6=$NET_FILE a7=$ACTION < $RPATH/L_final_plot.R

########################### 18. Calculate modules eigengene - SNP correlation   ###########################
RESULTS_TRANS=$EQTLRESULTS/3_output_matrixEQTL/results_trans.txt
RESULTS_CIS=$EQTLRESULTS/3_output_matrixEQTL/results_cis.txt
NETWORK_EIGENGENES=$EQTLRESULTS/4_downstream_matrixEQTL/co-expression_network/network_mods.txt
SNP_FILE=$EQTLRESULTS/snp_data_filtered.txt
OUTPUT_FOLDER=$EQTLRESULTS/4_downstream_matrixEQTL/co-expression_network_analysis
FIRST_SNP=20

R --no-save --args a1=$RESULTS_TRANS a2=$RESULTS_CIS a3=$NETWORK_EIGENGENES a4=$SNP_FILE a5=$OUTPUT_FOLDER a6=$FIRST_SNP  < $RPATH/Q_eigen-mostSNP_cor.R

########################### 19. eGenes vs non eGenes analysis: H2 and connectivity ###########################
CONN_FILE=$EQTLRESULTS/4_downstream_matrixEQTL/co-expression_network/gene_conn.txt
H2_FILE=$EQTLRESULTS/herit/H2.txt
OUTPUT_FILE=$EQTLRESULTS/herit/boxplot_H2+connectivity.pdf

RESULTS_TRANS=$EQTLRESULTS/3_output_matrixEQTL/results_trans.txt
RESULTS_CIS=$EQTLRESULTS/3_output_matrixEQTL/results_cis.txt
ALL_GENES=$EQTLRESULTS/expN_data_merged_filtered.txt
R --no-save --args a1=$CONN_FILE a2=$H2_FILE a3=$RESULTS_TRANS a4=$RESULTS_CIS a5=$ALL_GENES a6=$OUTPUT_FILE < $RPATH/R_h2+conn_analysis.R


############################################################################################################
#GENERAL LINEAR MODEL
############################################################################################################

########################### 20. General Linear Model ###########################
INPUT_FILE=$EQTLRESULTS/glm/glmTable.txt
INPUT_FILE_VG=$EQTLRESULTS/herit/var+H2.txt
OUTPUT_FOLDER=$EQTLRESULTS/glm/

R --no-save --args a1=$INPUT_FILE a2=$INPUT_FILE_VG a3=$OUTPUT_FOLDER < $RPATH/S_glm.R
