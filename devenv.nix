{ pkgs, ... }:

{
  languages.python = {
    enable = true;
    package = pkgs.python314;
    venv.enable = true;

    uv = {
      enable = true;
      sync.enable = true;
    };
  };

  env = {
    MARIMO_SKIP_UPDATE_CHECK = "1";
    PYTHONNOUSERSITE = "1";
  };
}
