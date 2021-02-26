#download mini conda from the following website
#https://docs.conda.io/en/latest/miniconda.html
#conda is installed to:
#~/miniconda3/condabin/conda

#download polyfun from the following website
#https://github.com/omerwe/polyfun
#polyfun is installed at /data/zhangh24/MR_MA/software/polyfun
#git clone https://github.com/omerwe/polyfun
#cd polyfun
#~/miniconda3/condabin/conda env create -f polyfun.yml
#source ~/miniconda3/etc/profile.d/conda.sh
#conda activate polyfun


#install susieR
#library(devtools)
install_github("stephenslab/susieR@0.8.0",build_vignettes=FALSE)
library(susieR)

#install LDscore 2.0
#wget http://www.christianbenner.com/ldstore_v2.0_x86_64.tgz
#tar xzvf ldstore_v2.0_x86_64.tgz








