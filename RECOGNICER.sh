#!/bin/sh

if [ $# -lt 3 ]; then
    echo ""
    echo 1>&2 Usage: $0 ["chip bedfile"] ["control bedfile"] ["FDR"]
    echo ""
    exit 1
fi

SRC_DIR=.
EX_DIR=$SRC_DIR/src

SAMPLEBED=$1
CONTROL=$2
FDR=$3

SAMPLE=Sample
SPECIES=hg18
OUT_DIR=.
WINDOW_SIZE=200
FRAGMENT_SIZE=150
STEP_SIZE=3
STEP_SCORE=2
GENOME_FRACTION=0.74
CHIPTHRESHOLD=1


echo ""
echo "###################################################"
echo "######            RECOGNICER v1.0            ######"
echo "###################################################"

echo ""
echo "Preprocessing the raw $SAMPLE file to remove redundant reads with threshold $CHIPTHRESHOLD..."
echo "python $EX_DIR/remove_redundant_reads.py -s $SPECIES -b $SAMPLEBED -t $CHIPTHRESHOLD -o $OUT_DIR/$SAMPLE-nonredundant.bed"
python $EX_DIR/remove_redundant_reads.py -s $SPECIES -b $SAMPLEBED -t $CHIPTHRESHOLD -o $OUT_DIR/$SAMPLE-nonredundant.bed

echo ""
echo "Generating bedgraph file..."
echo "python $EX_DIR/run-make-graph-file-by-chrom.py -s $SPECIES -b $OUT_DIR/$SAMPLE-nonredundant.bed -w $WINDOW_SIZE -i $FRAGMENT_SIZE -o $OUT_DIR/$SAMPLE-W$WINDOW_SIZE.bedgraph"
python $EX_DIR/run-make-graph-file-by-chrom.py -s $SPECIES -b $OUT_DIR/$SAMPLE-nonredundant.bed -w $WINDOW_SIZE -i $FRAGMENT_SIZE -o $OUT_DIR/$SAMPLE-W$WINDOW_SIZE.bedgraph

echo ""
echo "Coarse graining..."
echo "python $EX_DIR/coarsegraining-v9.py -s $SPECIES -b $OUT_DIR/$SAMPLE-W$WINDOW_SIZE.bedgraph -w $WINDOW_SIZE -g $STEP_SIZE -e $STEP_SCORE -t $GENOME_FRACTION -f $OUT_DIR/$SAMPLE.cgisland"
python $EX_DIR/coarsegraining-v9.py -s $SPECIES -b $OUT_DIR/$SAMPLE-W$WINDOW_SIZE.bedgraph -w $WINDOW_SIZE -g $STEP_SIZE -e $STEP_SCORE -t $GENOME_FRACTION -f $OUT_DIR/$SAMPLE.cgisland

echo ""
echo "Calculating significance of candidate islands comparing to the control..."
echo "python $EX_DIR/associate_tags_with_chip_and_control_w_fc_q_2013.py -s $SPECIES -f $FRAGMENT_SIZE -d $OUT_DIR/$SAMPLE.cgisland -b $CONTROL -a  $OUT_DIR/$SAMPLE-nonredundant.bed -o $OUT_DIR/$SAMPLE.cgsummary"
python $EX_DIR/associate_tags_with_chip_and_control_w_fc_q_2013.py -s $SPECIES -f $FRAGMENT_SIZE -d $OUT_DIR/$SAMPLE.cgisland -b $CONTROL -a  $OUT_DIR/$SAMPLE-nonredundant.bed -o $OUT_DIR/$SAMPLE.cgsummary

echo ""
echo "Finding significant islands using the FDR criterion..."
echo "python $EX_DIR/filter_islands_by_FDR_broadPeak.py -i $OUT_DIR/$SAMPLE.cgsummary -p $FDR -c 7 -o $OUT_DIR/${SAMPLE}-fdr${FDR}_broadPeak.bed"
python $EX_DIR/filter_islands_by_FDR_broadPeak.py -i $OUT_DIR/$SAMPLE.cgsummary -p $FDR -c 7 -o $OUT_DIR/${SAMPLE}-fdr${FDR}_broadPeak.bed

echo "Done!"