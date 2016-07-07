echo cd
cd

echo mkdir -p download
mkdir -p download

echo cd download
cd download

echo echo "Cached in $HOME/download :"
echo "Cached in $HOME/download :"

echo ls -l
ls -l


echo wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh
wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh

echo chmod +x miniconda.sh && ./miniconda.sh -b -p $HOME/miniconda
chmod +x miniconda.sh && ./miniconda.sh -b -p $HOME/miniconda

echo cd ..
cd ..

#echo export PATH="$HOME/miniconda/bin:$PATH"
#export PATH="$HOME/miniconda/bin:$PATH"


echo conda update --yes --quiet conda
conda update --yes --quiet conda

conda install python numpy scipy cython nose coverage matplotlib sphinx pillow

#echo conda create -n testenv --yes --quiet python numpy scipy cython nose coverage matplotlib sphinx pillow
#conda create -n testenv --yes --quiet python numpy scipy cython nose coverage matplotlib sphinx pillow

#echo source activate testenv
#source activate testenv

echo pip install .
pip install .
#sudo apt-get install python-numpy python-scipy
echo 
