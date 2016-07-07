cd
mkdir -p download
cd download
echo "Cached in $HOME/download :"
ls -l
wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh
chmod +x miniconda.sh && ./miniconda.sh -b -p $HOME/miniconda
cd ..
export PATH="$HOME/miniconda/bin:$PATH"
conda update --yes --quiet conda
conda create -n testenv --yes --quiet python numpy scipy cython nose coverage matplotlib sphinx pillow
source activate testenv     
#sudo apt-get install python-numpy python-scipy
pip install .
