{pkgs ? import <nixpkgs> {} }:

pkgs.mkShell {
    packages = [
        (pkgs.python3.withPackages(p: with p; [
            numpy
            ipympl
            ipywidgets
            matplotlib
            numpy
            pandas
            scipy
            sympy
        ]))
        pkgs.jupyter
        pkgs.nbstripout
    ];

}
