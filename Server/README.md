# GWAS_2
1. Симуляция гаплотипов через hapgen2 по данным 1000GP_Phase3

hapdir=~/1000GP_Phase3
for chr in `seq 1 22`; do
	gunzip $hapdir/1000GP_Phase3_chr${chr}.legend.gz
	gunzip $hapdir/1000GP_Phase3_chr${chr}.hap.gz

	dummyDL=`sed -n '2'p $hapdir/1000GP_Phase3_chr${chr}.legend | cut -d ' ' -f 2`
	~/hapgen2/hapgen2 -m $hapdir/genetic_map_chr${chr}_combined_b37.txt \
        -l $hapdir/1000GP_Phase3_chr${chr}.legend \
        -h $hapdir/1000GP_Phase3_chr${chr}.hap -o ~/hapgen_results/genotypes_chr${chr}_hapgen \
        -dl $dummyDL 0 0 0 -n 1000 0 -no_haps_output 

	rm $hapdir/1000GP_Phase3_chr${chr}.hap
done

2. Конвертация полученных гаплотипов в формат plink
for chr in `seq 1 22`; do
        plink --data genotypes_chr${chr}_hapgen.controls \
        --oxford-single-chr $chr \
        --make-bed \
        --out genotypes_chr${chr}_hapgen.controls
        
	echo -e "genotypes_chr${chr}_hapgen.controls" >> ~/hapgen_results/merge_total_list1
done

plink --merge-list file_list --make-bed --out genotypes_genome_hapgen.controls

3. Отбор вариант
~/plink/plink --bfile genotypes_genome_hapgen.controls --extract caspase_list.txt \
--make-bed --out genotypes_hapgen.controls.caspase

cut -f2 genotypes_genome_hapgen.controls.bim > snps.map 
shuf -n 985 snps.map > snps.subset.map
#руками добавил каузальные snp в общий датасет

~/plink/plink --bfile genotypes_genome_hapgen.controls --extract snps.subset.map \
--make-bed --out genotypes_subset_hapgen.controls

4. Построение фенотипа (см R_script2.R)
5. Ассоциация через plink
cat Y_caspase | awk 'BEGIN { FS="\t"; OFS="\t" } { $1=$1 "\t" $1 } 1' > Y_caspase1 
cat Y_caspase1 | awk '{FS="\t";OFS="\t"} {sub("id1_","id2_",$2)}1' > Y_caspase2
cat Y_caspase2 | tr -d '"' > Y_caspase3

~/GWAS/plink_linux_x86_64_20200219/plink --bfile ~/GWAS/Server/pathway/genotypes_subset_hapgen.controls \
--allow-no-sex --pheno ~/GWAS/Server/pathway/Y_caspase3 --all-pheno --assoc qt-means --out A_caspase


