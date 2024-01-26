This repository hosts the code corresponding to the research presented in the paper: [MELO et al. (2023) - Liquid–Liquid Equilibrium in Xylene Solubles (XS) Analysis of Polypropylene](https://onlinelibrary.wiley.com/doi/abs/10.1002/mren.202300029). The study employs a multicomponent Flory-Huggins model to elucidate the liquid–liquid equilibrium phenomenon within mixtures of polypropylene and xylene, focusing on the well-established xylene solubles (XS) test.

## Repository Structure

The majority of the code in this repository is written in Python:

- The `PolyLLE` class encompasses generic methods for:
  - Computing phase equilibrium compositions.
  - Estimating the Flory-Huggins parameter based on experimental data.
  - Plotting both experimental and calculated molecular weight distributions.
- The `Tompa` class, inheriting from `PolyLLE`, specifically implements the thermodynamic model described by [HEIDEMANN et al. (2006)](https://doi.org/10.1016/j.fluid.2005.12.030) and originally proposed by [TOMPA (1950)](https://pubs.rsc.org/en/content/articlelanding/1950/tf/tf9504600970).
- The `Kamide` class, also inheriting from `PolyLLE`, specifically implements the thermodynamic model outlined by [KAMIDE and DOBASHI (2000)](https://www.sciencedirect.com/book/9780444894304/physical-chemistry-of-polymer-solutions).
- The Jupyter notebooks generate results as presented in [MELO et al. (2023)](https://onlinelibrary.wiley.com/doi/abs/10.1002/mren.202300029). These notebooks effectively utilize the aforementioned classes to perform a test flash calculation and compute Flory-Huggins parameters from actual experimental data.

Moreover, the `Matlab` directory houses legacy code from our original implementation of the solution to this problem, as detailed in [MELO (2014)](http://hdl.handle.net/11422/7650) (my MSc dissertation in Portuguese).

## Citing

If this code has helped you in your research, consider citing:

```
@article{melo_liquidliquid_2023,
	title = {Liquid–{Liquid} {Equilibrium} in {Xylene} {Solubles} ({XS}) {Analysis} of {Polypropylene}},
	volume = {17},
	issn = {1862-832X, 1862-8338},
	url = {https://onlinelibrary.wiley.com/doi/10.1002/mren.202300029},
	doi = {10.1002/mren.202300029},
	number = {6},
	journal = {Macromolecular Reaction Engineering},
	author = {Melo, Afrânio and Pessoa, Fernando L. P. and Pinto, José Carlos},
	month = dec,
	year = {2023},
	pages = {2300029},
}
```

