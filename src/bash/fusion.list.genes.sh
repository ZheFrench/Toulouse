cd /home/jp/eclipse-workspace/Toulouse/data/results/genes_3_2/


cat pc9.genes.pattern1A.txt  pc9.genes.pattern2A.txt > pc9.genes.1A.2A.txt
cat pc9.genes.pattern1B.txt  pc9.genes.pattern2B.txt > pc9.genes.1B.2B.txt


cat h4.genes.pattern1A.txt  h4.genes.pattern2A.txt > h4.genes.1A.2A.txt
cat h4.genes.pattern1B.txt  h4.genes.pattern2B.txt > h4.genes.1B.2B.txt

wc -l pc9.genes.1A.2A.txt
wc -l pc9.genes.1B.2B.txt
wc -l h4.genes.1A.2A.txt
wc -l h4.genes.1B.2B.txt

#393
wc -l /home/jp/eclipse-workspace/Toulouse/data/gene.symbol.rnabindingprotein_human.csv 
grep -f /home/jp/eclipse-workspace/Toulouse/data/gene.symbol.rnabindingprotein_human.csv -Hrn /home/jp/eclipse-workspace/Toulouse/data/results/genes_3_2/ | grep genes.pattern | sed -r 's/.genes.pattern/\t/' | sed -r s'/.txt:[0-9]+:/\t/' | sed -r s'/.*\///' > /home/jp/eclipse-workspace/Toulouse/data/results/majorRegulators/rna.tsv
awk -F "\t" '{print $3}' /home/jp/eclipse-workspace/Toulouse/data/results/majorRegulators/rna.tsv > ~/Desktop/eclipse-workspace/Toulouse/data/rna.txt

#1640
wc -l /home/jp/eclipse-workspace/Toulouse/data/gene.symbol.dnabindingprotein_human.csv 
grep -f /home/jp/eclipse-workspace/Toulouse/data/gene.symbol.dnabindingprotein_human.csv -Hrn /home/jp/eclipse-workspace/Toulouse/data/results/genes_3_2/ | grep genes.pattern | sed -r 's/.genes.pattern/\t/' | sed -r s'/.txt:[0-9]+:/\t/' | sed -r s'/.*\///' > /home/jp/eclipse-workspace/Toulouse/data/results/majorRegulators/dna.tsv
awk -F "\t" '{print $3}' /home/jp/eclipse-workspace/Toulouse/data/results/majorRegulators/dna.tsv > ~/Desktop/eclipse-workspace/Toulouse/data/dna.txt

#430
wc -l /home/jp/eclipse-workspace/Toulouse/data/gene.symbol.chromatinRemodelers_human.csv
grep -f /home/jp/eclipse-workspace/Toulouse/data/gene.symbol.chromatinRemodelers_human.csv -Hrn /home/jp/eclipse-workspace/Toulouse/data/results/genes_3_2/ | grep genes.pattern | sed -r 's/.genes.pattern/\t/' | sed -r s'/.txt:[0-9]+:/\t/' | sed -r s'/.*\///' > /home/jp/eclipse-workspace/Toulouse/data/results/majorRegulators/chrom.tsv
awk -F "\t" '{print $3}' /home/jp/eclipse-workspace/Toulouse/data/results/majorRegulators/chrom.tsv > ~/Desktop/eclipse-workspace/Toulouse/data/chrom.txt
