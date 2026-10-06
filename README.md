# Meet-in-the-MidL

This repository hosts the codebase accompanying the article *Meet-in-the-midL: the complexity of graph walking on compositions of isogeny graphs*.

### Requirements
A recent version of Sagemath (>= 10.5) is required. The computations are accelerated using pari.gp, which is bundled with Sagemath by default.

### Configuration
By default, the experiments utilise the library [Isogeny](https://github.com/JamesRickards-Canada/Isogeny/) for computing the supersingular $L$-isogeny graph(s) $G_{p, L}$. Full guidance for configuring this package may be found within this directory.

In short, run
```bash
cd Isogeny
./configure
```
This will prompt for a path to the file `pari.cfg` bundled with a pari installation. It is recommended to point this to the pari installation bundled with Sagemath.
```bash
make
cd ..
ln -s Isogeny/isogeny.gp
ln -s Isogeny/libisogeny-X-Y-Z.so
ln -s Isogeny/modpolys/
```
where `X, Y, Z` should be substituted to match the name of the file created by `make`.

### Usage
The script `main.py` generates an example plot of a cumulative distribution.