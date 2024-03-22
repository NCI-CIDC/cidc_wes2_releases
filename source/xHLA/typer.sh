#!/bin/bash
# Author: Xie Chao
set -eu -o pipefail

[[ $# -lt 2 ]] && {
    echo "usage: $(basename "$0") S3://path.bam sample_id [delete]";
    exit 1;
}

##THIS line specifies the directory for locating the xHLA scripts
BIN="`dirname \"$0\"`"
##echo "BIN VALUE:",$BIN

S3=$1
ID=$2
HLA_REF=$3
#OUT needs to be an absolute path
OUT=xhla-$ID #old command as relative path
OUTPATH=$4
DELETE=false
FULL=false

echo $OUT

while test $# -gt 0
do
    case "$1" in
        delete) DELETE=true
        ;;
        full) FULL=true
        ;;
    esac
    shift
done

echo "typer.sh parameters: DELETE=$DELETE FULL=$FULL"

mkdir -p $OUTPATH/$OUT
TEMP=temp-$RANDOM-$RANDOM-$RANDOM

echo "Extracting reads from S3"

samtools view -u $S3 chr6:29886751-33090696 | samtools view -L $HLA_REF/hla.bed - > ${TEMP}.sam
$BIN/preprocess.pl ${TEMP}.sam | gzip > $OUTPATH/$OUT/$ID.fq.gz
rm ${TEMP}.sam
echo "Aligning "$OUTPATH/$OUT/$ID.fq.gz " reads to IMGT database"

if [ "$FULL" = true ]; then
    $BIN/align.pl $OUTPATH/$OUT/${ID}.fq.gz $OUTPATH/$OUT/${ID}.tsv $HLA_REF full
else
    $BIN/align.pl $OUTPATH/$OUT/${ID}.fq.gz $OUTPATH/$OUT/${ID}.tsv $HLA_REF
fi

echo "Typing"
echo "running typing command: " $BIN/typing.r $OUTPATH/$OUT/${ID}.tsv $OUTPATH/$OUT/${ID}.hla.tsv $HLA_REF
$BIN/typing.r $OUTPATH/$OUT/${ID}.tsv $OUTPATH/$OUT/${ID}.hla.tsv $HLA_REF

echo "Reporting"
echo "running reporting command:" $BIN/report.py -in $OUTPATH/$OUT/${ID}.hla.tsv -out $OUTPATH/$OUT/${ID}.json -subject $ID -sample $ID
$BIN/report.py -in $OUTPATH/$OUT/${ID}.hla.tsv -out $OUTPATH/$OUT/${ID}.json -subject $ID -sample $ID

if [ "$FULL" = true ]; then
    echo "running command:" $BIN/full.r $OUTPATH/$OUT/${ID}.tsv.dna $OUTPATH/$OUT/${ID}.hla.tsv $OUTPATH/$OUT/${ID}.hla.full
    $BIN/full.r $OUTPATH/$OUT/${ID}.tsv.dna $OUTPATH/$OUT/${ID}.hla.tsv $OUTPATH/$OUT/${ID}.hla.full
fi

# Clean up
if [ "$DELETE" = true ]
then
	rm $OUTPATH/$OUT/${ID}.tsv
	rm $OUTPATH/$OUT/${ID}.fq.gz
	rm $OUTPATH/$OUT/${ID}.hla
fi
