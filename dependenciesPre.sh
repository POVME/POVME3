echo rm -rf virtualenvs venv .pyenv
rm -rf virtualenvs venv .pyenv
echo
echo

echo deactivate
deactivate
echo
echo

echo cd
cd
echo
echo

echo mkdir -p download
mkdir -p download
echo
echo

echo cd download
cd download
echo
echo

echo echo "Cached in $HOME/download :"
echo "Cached in $HOME/download :"
echo
echo

echo ls -l
ls -l
echo
echo


echo wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh
wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh
echo
echo

echo chmod +x miniconda.sh && ./miniconda.sh -b -p $HOME/miniconda
chmod +x miniconda.sh && ./miniconda.sh -b -p $HOME/miniconda
echo
echo

echo cd ..
cd ..
echo
echo

#echo export PATH="$HOME/miniconda/bin:$PATH"
#export PATH="$HOME/miniconda/bin:$PATH"


echo conda update --yes --quiet conda
conda update --yes --quiet conda
echo
echo

echo conda install --yes python numpy scipy 
conda install --yes python numpy scipy 
echo
echo

#echo conda create -n testenv --yes --quiet python numpy scipy cython nose coverage matplotlib sphinx pillow
#conda create -n testenv --yes --quiet python numpy scipy cython nose coverage matplotlib sphinx pillow

#echo source activate testenv
#source activate testenv

#echo pip install .
#pip install .
#echo
#echo
#sudo apt-get install python-numpy python-scipy

echo cd
cd
echo
echo

echo cd POVME
cd POVME
echo
echo

echo python setup.py install
python setup.py install
echo 
echo
