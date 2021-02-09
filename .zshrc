export PATH=/usr/bin/python:$PATH

#export FITPACK=/Users/patrickbarry/work/JAM/fitpack
#export FITPACK=/Users/patrickbarry/work/JAM/fitpack-resum
export FITPACK=/Users/patrickbarry/work/JAM/fitpack2
export ML4JAM=/Users/patrickbarry/work/JAM/ml4jam

export PATH=$FITPACK:$PATH
export PYTHONPATH=$FITPACK:$PYTHONPATH

#--for python2
export PATH=/usr/local/Cellar/lhapdf/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/Cellar/lhapdf/lib:$LD_LIBRARY_PATH
export PYTHONPATH=/usr/local/Cellar/lhapdf/lib/python2.7/site-packages/:$PYTHONPATH

export LHAPDF_DATA_PATH=/usr/local/Cellar/lhapdf/share/LHAPDF:$LHAPDF_DATA_PATH

#--for python3
export PATH=/usr/local/Cellar/lhapdf3/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/Cellar/lhapdf3/lib:$LD_LIBRARY_PATH
export PYTHONPATH=/usr/local/Cellar/lhapdf3/lib/python3.7/site-packages/:$PYTHONPATH

export LHAPDF_DATA_PATH=/usr/local/Cellar/lhapdf3/share/LHAPDF:$LHAPDF_DATA_PATH

# export PATH="/Users/patrickbarry/opt/anaconda2/bin:$PATH"  # commented out by conda initialize

#--to see if things work for lhapdf
#export PATH=$HOME/lhapdf/bin:$PATH
#export LD_LIBRARY_PATH=$HOME/lhapdf/lib:$LD_LIBRARY_PATH
#export PYTHONPATH=$HOME/lhapdf/lib/python2.7/site-packages/:$PYTHONPATH
#export LHAPDF_DATA_PATH=$HOME/lhapdf/share/LHAPDF:$LHAPDF_DATA_PATH


export CLICOLOR=1

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/Users/patrickbarry/opt/anaconda2/bin/conda' 'shell.zsh' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/Users/patrickbarry/opt/anaconda2/etc/profile.d/conda.sh" ]; then
        . "/Users/patrickbarry/opt/anaconda2/etc/profile.d/conda.sh"
    else
        export PATH="/Users/patrickbarry/opt/anaconda2/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda deactivate
