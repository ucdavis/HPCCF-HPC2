#%Module1.0

set name "Miniconda3"
set dir "conda3"
set version "4.X"
set os "ucdhpc-20.04"

set url "https://docs.conda.io/en/latest/miniconda.html"
set tags "programming language"

set note "Miniconda is a free minimal installer for conda. It is a small, bootstrap version of Anaconda that includes only conda, Python, the packages they depend on, and a small number of other useful packages, including pip, zlib and a few others. Use the conda install command to install 720+ additional conda packages from the Anaconda repository."

module-whatis "This module loads $name version $version. $note"

if { [info exists env(CONDA_EXE)] && [module-info command load] } {
    puts stderr "ERROR: CONDA_EXE is currently defined: $env(CONDA_EXE)."
    puts stderr "This module will almost certainly interfere with your conda installation."
    puts stderr "Remove existing conda installation and its shell hooks from your PATH before proceeding."
    exit
}

if { [info exists env(CONDA_PREFIX)] && [module-info command unload] } {
    puts stderr "ERROR: CONDA_PREFIX is currently defined: $env(CONDA_PREFIX)."
    puts stderr "This means you have a conda environment active."
    puts stderr "Make sure to `conda deactivate` before you unload this module."
    exit
}

set swroot "/software"
prepend-path PATH            $swroot/$dir/$version/bin
prepend-path PATH            $swroot/$dir/$version/condabin
prepend-path MANPATH         $swroot/$dir/$version/share/man

# https://setuptools.pypa.io/en/latest/deprecated/distutils-legacy.html
setenv SETUPTOOLS_USE_DISTUTILS stdlib
setenv CONDA_MOD_ROOT $swroot/$dir/$version/

switch -- [module-info shell] {
    bash {
        source-sh bash $env(CONDA_MOD_ROOT)/shell-setup.bash
    }
    sh {
        source-sh bash $env(CONDA_MOD_ROOT)/shell-setup.bash
    }
    zsh {
        source-sh zsh $env(CONDA_MOD_ROOT)/shell-setup.zsh
    }
    csh {
        source-sh tcsh $env(CONDA_MOD_ROOT)/shell-setup.tcsh
    }
}

if { [ module-info mode load ] } {
    puts stderr "Module $dir/$version loaded."
}
