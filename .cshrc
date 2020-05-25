

#ident  "@(#).cshrc     ver 1.0     Aug 20, 1996"
# Default user .cshrc file.
#
# This file is executed each time a shell is started.
# This includes the execution of shell scripts.


#####
# Source the site-wide syscshrc file.
# The syscshrc file defines some needed aliases (setup amd unsetup)
# and environment variables (PATH and MANPATH).  This line
# should not be deleted.  You do, however, have a choice of
# syscshrc files.  Uncomment the one that you prefer.
#####
source /site/env/syscshrc       # Searches /usr/local/bin first.
#source /site/env/syscshrc.alt   # Searches /usr/local/bin last.


#####
# Set up the shell environment.  You may comment/uncomment
# the following entries to meet your needs.
#####
# Number of commands to save in history list.
set history=50
 
# Number of commands to save in ~/.history upon logout.
set savehist=50

# Notify user of completed jobs right away, instead of waiting
# for the next prompt.
#set notify

# Don't redirect output to an existing file.
# CAD NOTE!  This must be commented out for proper ME10 functionality!!
set noclobber

# Set the file creation mode mask (default permissions for newly created files).
umask 022


#####
# Define your aliases.
#####
#alias       h       history
#alias       d       dirs
#alias       pd      pushd
#alias       pd2     pushd +2
#alias       po      popd
#alias       m       more
#alias       ls      'ls -F'


#####
# Define your default printer.
#####

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
if ( -f "/work/JAM/apps/anaconda3/etc/profile.d/conda.csh" ) then
    source "/work/JAM/apps/anaconda3/etc/profile.d/conda.csh"
else
    setenv PATH="/work/JAM/apps/anaconda3/bin:$PATH"
endif
# <<< conda initialize <<<

## >>> conda initialize >>>
## !! Contents within this block are managed by 'conda init' !!
#if ( -f "/work/JAM/apps/anaconda2/etc/profile.d/conda.csh" ) then
#    source "/work/JAM/apps/anaconda2/etc/profile.d/conda.csh"
#else
#    setenv PATH="/work/JAM/apps/anaconda2/bin:$PATH"
#endif
## <<< conda initialize <<<
#
#####
# User specific additions should be added below.
#####
setenv TERM xterm-color
#setenv TERM xterm
setenv PATH ~/apps/bin:${PATH}
setenv PATH /site/bin:${PATH}
setenv PATH /apps/tmux/bin:${PATH}
setenv PATH /apps/vim/bin:${PATH}
unsetenv SSH_ASKPASS
alias desk ssh thypc21
alias iml ssh iml19g01
alias wgetncc wget --no-check-certificate
setenv LSCOLORS ExFxCxDxBxegedabagacad
alias c "clear"
alias ls 'ls --color=always --ignore="*.pyc"  --ignore="*__.py"  --ignore="__pycache__"  '
#alias solid 'cd /lustre/expphy/work/halla/solid/nsato'
#alias cteqx 'cd /u/group/cteqX'
#alias eva 'cd /lustre/expphy/work/hallb/clas12/avakian/eva'
#alias jam 'cd /lustre/expphy/volatile/theory/cteqX'
#alias web 'cd /u/site/www/html/theory'
#alias jamdold 'cd /lustre/volatile/JAM'
alias jamd 'cd /w/general-scifs17exp/JAM'
#source ~/.tmux-completion.tcsh 
alias ml 'ssh iml19g01'


#--anaconda
setenv PATH /work/JAM/apps/anaconda3/bin:${PATH}
#setenv PATH /work/JAM/apps/anaconda2/bin:${PATH}
#setenv PATH /u/apps/anaconda/python2/anaconda201812/bin:${PATH}
#setenv PATH /u/apps/anaconda/python3/anaconda201812/bin:${PATH}

#--pypy
setenv PATH /enp/w/general-scifs17exp/JAM/apps/pypy-7.1.1-linux_x86_64-portable/bin:${PATH}

#--hack for CXXABI_1.3.8
#setenv LD_LIBRARY_PATH /work/JAM/apps/lib/:${LD_LIBRARY_PATH}

#--jam1d 
setenv FITPACK /work/JAM/barryp/JAM/fitpack/
setenv ML4JAM /work/JAM/barryp/JAM/ml4jam
#setenv PATH ${FITPACK}/bin:${PATH}
setenv PATH ${FITPACK}:${PATH}
setenv PYTHONPATH ${FITPACK}

#setenv STFITPACK /work/JAM/nsato/qcdhub/stf-fitter/
#setenv PATH ${STFITPACK}/bin:${PATH}
 
#--lhadpf
#setenv PATH /work/JAM/apps/lhapdf2/bin:${PATH}
#setenv PYTHONPATH /work/JAM/apps/lhapdf2/lib/python2.7/site-packages/
#setenv LD_LIBRARY_PATH /work/JAM/apps/lhapdf2/lib#:${LD_LIBRARY_PATH}
setenv PATH /work/JAM/apps/lhapdf2/bin:${PATH}
setenv LD_LIBRARY_PATH /work/JAM/apps/lhapdf2/lib
setenv PYTHONPATH /work/JAM/apps/lhapdf2/lib/python2.7/site-packages/:${PYTHONPATH}
setenv LHAPDF_DATA_PATH /work/JAM/apps/lhapdf2/share/LHAPDF/


#--HIJET 
#setenv HIJPACK /work/JAM/nsato/qcdhub/hijetpack/
#setenv PATH ${HIJPACK}/bin:${PATH}
#setenv PYTHONPATH ${HIJPACK}:${PYTHONPATH}


#--slurm
#setenv PATH /site/scicomp/auger-slurm/bin:${PATH}

#--slurm
#setenv PATH /work/JAM/nsato/myapps:${PATH}

#--singularity
#setenv PATH /apps/singularity/bin/:${PATH}

#--fastjet
#setenv PATH /work/JAM/apps/fastjet/bin/:${PATH}


#module load gcc_7.2.0

#setenv  HDF5_USE_FILE_LOCKING FALSE



#alias 12s1 ssh qcd12s0241
#alias 12s2 ssh qcd12s0242







