# Monte-Carlo for NIR-II Imaging

Monte-Carlo simulation for NIR-II Imaging.

This project uses the method proposed by [Wang et al.](https://doi.org/10.1016/0169-2607(95)01640-F) to simulating photon propagation with the Monte Carlo model. The number of photons is set to 100,000 for each wave length.

## Usage

Python 3 is required.

Install [pipx](https://github.com/pypa/pipx) >=1.4.2:

```sh
python -m pip install -U pipx
```

Run the Monte-Carlo simulation:

```sh
pipx run run.py > output.txt
```

The output files `output.txt`, `pDW`, and `pz` has been given in this repository.

To process your own data, modify `parameters.csv` file.

## Notice

Derived from [https://github.com/titonmoy/Monte-Carlo](https://github.com/titonmoy/Monte-Carlo) (authored by Thamidul Islam Tonmoy) under [MIT license](LICENSE).

## License

[MIT license](LICENSE)
