{ config, pkgs, ... }:

{
  packages = [
    pkgs.blast
    pkgs.mafft
  ];

  languages.python = {
    enable = true;
    package = pkgs.python314;
    venv.enable = true;

    uv = {
      enable = true;
      sync.enable = true;
    };
  };

  tasks."project:sync-conda-environment" = {
    description = "Synchronize the generated Conda environment specification";
    exec = "python scripts/sync_conda_environment.py";
    cwd = config.git.root;
    after = [ "devenv:python:uv" ];
    before = [ "devenv:enterShell" ];
    execIfModified = [
      "pyproject.toml"
      "uv.lock"
      "scripts/sync_conda_environment.py"
      "environment.yml"
    ];
  };

  env = {
    MARIMO_SKIP_UPDATE_CHECK = "1";
    PYTHONNOUSERSITE = "1";
  };
}
