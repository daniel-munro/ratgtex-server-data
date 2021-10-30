set -euxo pipefail

mkdir -p $2
mkdir -p $2/eqtl

echo "Assembling eQTL files"
Rscript eqtl_files.R $@

echo "Making tissue info table"
Rscript tissueInfo.R $@
