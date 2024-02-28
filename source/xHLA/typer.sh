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
OUT=xhla-$ID
DELETE=false
FULL=false

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

mkdir -p $OUT
TEMP=temp-$RANDOM-$RANDOM-$RANDOM

echo "Extracting reads from S3"

samtools view -u $S3 chr6:29886751-33090696 | samtools view -L $HLA_REF/hla.bed - > ${TEMP}.sam
$BIN/preprocess.pl ${TEMP}.sam | gzip > $OUT/$ID.fq.gz
rm ${TEMP}.sam
echo "Aligning reads to IMGT database"
if [ "$FULL" = true ]; then
    $BIN/align.pl $PWD/$OUT/${ID}.fq.gz $OUT/${ID}.tsv $HLA_REF full
else
    $BIN/align.pl $PWD/$OUT/${ID}.fq.gz $OUT/${ID}.tsv $HLA_REF
fi

echo "Typing"
echo "running command: " $BIN/typing.r $OUT/${ID}.tsv $OUT/${ID}.hla
$BIN/typing.r $OUT/${ID}.tsv $OUT/${ID}.hla

echo "Reporting"
$BIN/report.py -in $OUT/${ID}.hla -out $OUT/${ID}.json -subject $ID -sample $ID

if [ "$FULL" = true ]; then
    $BIN/full.r $OUT/${ID}.tsv.dna $OUT/${ID}.hla $OUT/${ID}.hla.full
fi

# Clean up
if [ "$DELETE" = true ]
then
	rm $OUT/${ID}.tsv
	rm $OUT/${ID}.fq.gz
	rm $OUT/${ID}.hla
fi
