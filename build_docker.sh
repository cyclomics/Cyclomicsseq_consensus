TMP_NAME=0.1.1

docker build -t cyclomics/cyclomicsseq_consensus:$TMP_NAME .
echo "Container created: cyclomics/cyclomicsseq_consensus:$TMP_NAME"
docker push cyclomics/cyclomicsseq_consensus:$TMP_NAME
echo "Container pushed: cyclomics/cyclomicsseq_consensus:$TMP_NAME"