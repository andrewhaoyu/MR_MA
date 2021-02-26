
for chr in {1..23}
do 
/data/zhangh24/software/plink2 --bfile /data/zhangh24/KG.plink/EUR/chr_all --chr $chr --make-bed --out /data/zhangh24/KG.plink/EUR/chr_${chr}; \
done