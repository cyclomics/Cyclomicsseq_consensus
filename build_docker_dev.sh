TMP_NAME=$(date +'%d%m%Y')

docker build -t cyclomics/cyclomicsseq_consensus:dev-$TMP_NAME .
echo "Container created cyclomics/cyclomicsseq_consensus:dev-$TMP_NAME"