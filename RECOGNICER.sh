#!/bin/sh

if [ $# -lt 3 ]; then
    echo ""
    echo 1>&2 Usage: $0 ["bed file"] ["control file"] ["p-value"] 
    echo ""
    exit 1
fi

WINDOW_SIZE=200
FRAGMENT_SIZE=150
STEP_SIZE=3
STEP_SCORE=2

SAMPLEBED=$1
SAMPLE=${SAMPLEBED%.*}

CONTROL=$2
PVALUE=$3

SPECIES=hg18

SAMPLE_DIR=../CD4_H3K27me3
CONTROL_DIR=../CD4_H3K27me3
OUT_DIR=.
EX_DIR=/cluster1/czang/coarse-graining-code/cg/Modules
 
echo "python $EX_DIR/run-make-graph-file-by-chrom.py -s $SPECIES -b $SAMPLE_DIR/$SAMPLEBED -w $WINDOW_SIZE -i $FRAGMENT_SIZE -o $OUT_DIR/$SAMPLE-W$WINDOW_SIZE.graph"
python $EX_DIR/run-make-graph-file-by-chrom.py -s $SPECIES -b $SAMPLE_DIR/$SAMPLEBED -w $WINDOW_SIZE -i $FRAGMENT_SIZE -o $OUT_DIR/$SAMPLE-W$WINDOW_SIZE.graph

echo "python $EX_DIR/coarsegraining-v9.py -s $SPECIES -b $OUT_DIR/$SAMPLE-W$WINDOW_SIZE.graph -w $WINDOW_SIZE -g $STEP_SIZE -e $STEP_SCORE -t $GENOME_FRACTION -f $OUT_DIR/$SAMPLE.cgisland"
python $EX_DIR/coarsegraining-v9.py -s $SPECIES -b $OUT_DIR/$SAMPLE-W$WINDOW_SIZE.graph -w $WINDOW_SIZE -g $STEP_SIZE -e $STEP_SCORE -t $GENOME_FRACTION -f $OUT_DIR/$SAMPLE.cgisland

echo "python $EX_DIR/associate_tags_with_chip_and_control_w_fc_q_2013.py -s $SPECIES -f $FRAGMENT_SIZE -d $OUT_DIR/$SAMPLE.cgisland -b $CONTROL_DIR/$CONTROL -a  $SAMPLE_DIR/$SAMPLEBED -o $OUT_DIR/$SAMPLE.cgsummary3"
python $EX_DIR/associate_tags_with_chip_and_control_w_fc_q_2013.py -s $SPECIES -f $FRAGMENT_SIZE -d $OUT_DIR/$SAMPLE.cgisland -b $CONTROL_DIR/$CONTROL -a  $SAMPLE_DIR/$SAMPLEBED -o $OUT_DIR/$SAMPLE.cgsummary

echo "python $EX_DIR/find_significant_islands.py -i $OUT_DIR/$SAMPLE.cgsummary3 -p $PVALUE -o $OUT_DIR/${SAMPLE}-P$PVALUE.bed"
python $EX_DIR/find_significant_islands.py -i $OUT_DIR/$SAMPLE.cgsummary -p $PVALUE -o $OUT_DIR/${SAMPLE}-P$PVALUE.bed

