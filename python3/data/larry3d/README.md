# Data for the Larry 3-D module

Travel time data sets and reference models used by the 3-D body wave
tomography module. Each subdirectory was kindly provided by the
authors named in its notes; please cite the papers listed there if
you use the data.

- `hrvdata/`: P wave summary (binned) delay times from Antolik,
  Ekström and Dziewonski (2001, *Pure Appl. Geophys.*), see
  `read.me` and the included paper. Default data set of the module.
- `houser_s/`: S wave data of Houser et al. (2008, *Geophys. J. Int.*),
  see the included paper.
- `ritsema_s/`: S wave measurements of Ritsema et al. (2011,
  *Geophys. J. Int.*), see `readme.txt`.
- `aniprmc`: the anisotropic PREM reference model used for the
  inversions.

`data.*.bin` files are generated from the `*.ascii` files by
`data_a2b` during installation; `conv.log` and `error.log` are the
logs of that conversion.
