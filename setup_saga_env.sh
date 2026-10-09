#!/usr/bin/env bash

# Source this file from any directory:
#   source /cluster/projects/nn10058k/leo/gref4hsi/setup_saga_env.sh

_gref4hsi_repo="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_gref4hsi_env="/cluster/work/users/leoguenz/venvs/gref4hsi_env"

if [[ ! -x "${_gref4hsi_env}/bin/python" ]]; then
    echo "gref4hsi environment not found: ${_gref4hsi_env}" >&2
    unset _gref4hsi_repo _gref4hsi_env
    return 1
fi

# Avoid packages inherited from ~/.local or unrelated ROS/workstation setups.
unset PYTHONPATH
export PYTHONNOUSERSITE=1

export CONDA_PREFIX="${_gref4hsi_env}"
export CONDA_DEFAULT_ENV="${_gref4hsi_env}"
export PATH="${_gref4hsi_env}/bin:${PATH}"
export LD_LIBRARY_PATH="${_gref4hsi_env}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

export GREF4HSI_ROOT="${_gref4hsi_repo}"
export PIP_CACHE_DIR="/cluster/work/users/leoguenz/pip_cache"
export XDG_CACHE_HOME="/cluster/work/users/leoguenz/.cache"
export GDAL_DATA="${_gref4hsi_env}/share/gdal"
export PROJ_DATA="${_gref4hsi_env}/share/proj"

# Safe defaults for compute nodes without an X display. Override after sourcing
# when an interactive GUI is explicitly required.
export MPLBACKEND="${MPLBACKEND:-Agg}"
export PYVISTA_OFF_SCREEN="${PYVISTA_OFF_SCREEN:-true}"
export QT_QPA_PLATFORM="${QT_QPA_PLATFORM:-offscreen}"

cd "${GREF4HSI_ROOT}"
unset _gref4hsi_repo _gref4hsi_env
