# mlst_scheme_evaluation

Enter the reproducible development shell and open the notebook with:

```console
devenv shell
marimo edit schema_evaluation.py
```

Python, uv, marimo, and the notebook dependencies are declared in
`pyproject.toml` and locked in `uv.lock`. `devenv shell` runs `uv sync`
automatically and activates the devenv-managed virtual environment. Use
`uv add <package>` to add or update a dependency.

On a machine without a usable NVIDIA or AMD GPU, marimo 0.25.x may emit a
warning that it cannot collect GPU statistics. It is harmless. To suppress
that upstream warning when desired, start marimo with:

```console
marimo --log-level error edit schema_evaluation.py
```
