# GWAS
Порядок команд в пайплайне:

1.Отбираю данные по 3м хромосомам и сливаю в один файл

~/GWAS/plink_linux_x86_64_20200219/plink --bfile EUR_1KG_nodups --chr 5 --make-bed --out EUR_1KG_chr5_nodup

~/GWAS/plink_linux_x86_64_20200219/plink --bfile EUR_1KG_nodups --chr 7 --make-bed --out EUR_1KG_chr7_nodup

~/GWAS/plink_linux_x86_64_20200219/plink --bfile EUR_1KG_nodups --chr 11 --make-bed --out EUR_1KG_chr11_nodup

~/GWAS/plink_linux_x86_64_20200219/plink --merge-list pruned_list --make-bed --allow-no-sex --out EUR_1KG_genome_nodup

2.Строю фенотип в PhenotypeSimulator (см код в R_script_manual.R)

3.GWAS ассоциация

~/GWAS/plink_linux_x86_64_20200219/plink --bfile ~/GWAS/EUR_1KG/EUR_1KG_genome_nodup --allow-no-sex --pheno ~/GWAS/Pheno2/simulation/Ysim_plink2.txt --tail-pheno 0 --assoc --out plink3

#lsea_table - реформатированный аутпут plink3.assoc

4.Анализ LSEA

python3 ~/GWAS/ukb_phewas/LSEA/LSEA.py -af ~/GWAS/lsea_table -sn ~/GWAS/snpEff_latest_core/snpEff/scripts/snpEff -pld ~/GWAS/plink_linux_x86_64_20200219/plink -bf ~/GWAS/EUR_1KG/EUR_1KG_genome_nodup -p

#На LSEA пока остановился, так как проблемы с snpEff
