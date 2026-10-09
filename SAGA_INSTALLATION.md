# Saga environment

The persistent environment is installed at:

```text
/cluster/work/users/leoguenz/venvs/gref4hsi_env
```

Activate it from any shell or Slurm job with:

```bash
source /cluster/projects/nn10058k/leo/gref4hsi/setup_saga_env.sh
```

The helper activates the environment, isolates it from user-site and inherited
`PYTHONPATH` packages, sets the GDAL and PROJ data directories, configures
writable cache locations, selects headless-safe plotting defaults, and changes
to the repository root. The repository is installed editable, so source edits
are immediately visible without reinstalling.

For an interactive graphical session, override the headless defaults after
activation as needed:

```bash
unset PYVISTA_OFF_SCREEN
unset QT_QPA_PLATFORM
unset MPLBACKEND
```

`requirements.txt` also lists the unmaintained `pyembree` package. It is not
installed because its build is incompatible with the current packaging stack
and it conflicts with current `trimesh`. The maintained `embreex` replacement
is installed and provides accelerated ray tracing.
