set -e
cd /home/dmin/Grants/Nitrogen_Cycle/ncycle-pipeline/comparators/nmag
rm -rf dram_out dram_distill
DRAM.py annotate_genes -i 'proteomes/*.faa' -o dram_out --threads 8
DRAM.py distill -i dram_out/annotations.tsv -o dram_distill
echo "DRAM_DONE"
