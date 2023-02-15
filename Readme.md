# APCMA modulation/demodulation block for GNU Radio and USRP


## how to build
OS : ubutnu20.04
```shell
git clone https://github.com/Atsushi0803/gr-apcma_sdr.git
cd ./gr-apcma_sdr
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash ./Miniconda3-latest-Linux-x86_64.sh
conda env create -f ./environment.xml
mkdir build
cd build
cmake ../
sudo make install
```
