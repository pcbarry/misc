# unpolarized Jet data from Tevatron and RHIC

## data tables

| index | ref                    | normalization | collision | year                                   | collaboration    | status |
| ----- | -----                  | -----         | -----     | -----                                  | -----            | -----  |
| 10001 | [HEP Data][link.10001] | no            | `ppb`     | 2004 to 2005 run II data               | D0               | ready  |
| 10002 | [HEP Data][link.10002] | yes           | `ppb`     | 2007 paper on 2002 to 2006 run II data | CDF              | ready  |
| 10003 | [drupal][link.10003]   | no            | `pp`      | 2006 paper on 2003 to 2004 data        | STAR MB          | ready  |
| 10004 | [drupal][link.10004]   | no            | `pp`      | 2006 paper on 2003 to 2004 data        | STAR HT          | ready  |

<br/>

1. dataset [10099][link.10099] is the cone algorithm version of dataset 10002

2. dataset 10082 is a copy of dataset 10099 except

    - it treats systematic uncertainties that are only correlated within single $\left| \eta \right|$ bin as fully correlated over all bins

3. dataset 10092 is a copy of dataset 10099 except

    - it adds all systematic uncertainties in 10099 in quadrature to give one single overall systematic uncertainty
    - it treats this one single systematic uncertainty as uncorrelated

3. dataset 10072 is a copy of dataset 10002 except

    - it separates the total systematic uncertainty into components
    - it treats the systematic uncertainties correlated as described in the paper

4. dataset 19999 is a test dataset

[link.10001]: https://www.hepdata.net/record/ins779574
[link.10002]: https://www.hepdata.net/record/ins743342
[link.10003]: https://drupal.star.bnl.gov/STAR/files/starpublications/68/data.html
[link.10004]: https://drupal.star.bnl.gov/STAR/files/starpublications/68/data.html
[link.10099]: https://www.hepdata.net/record/ins790693

## observables

- $\frac{\mathrm{d}^2 \sigma}{\mathrm{d} p_T \mathrm{d} y}$
- $\frac{\mathrm{d}^2 \sigma}{2 \pi \mathrm{d} p_T \mathrm{d} y}$

## headers

- `idx`: indices
- `col`: collaboration
- `particles-in`: `pp` for proton proton collision and `ppb` for proton anti proton collision
- `RS`: $\sqrt{s}$ in GeV
- `pt-min`: minimum $p_T$ in GeV
- `pt-max`: maximum $p_T$ in GeV
- `pT`: average $p_T$ in GeV
- `tau`: $\frac{2 p_t}{\sqrt{s}}$
- `eta-abs-min`: minimum $\left| \eta \right|$
- `eta-abs-max`: maximum $\left| \eta \right|$
- `eta-min`: minimum $\eta$
- `eta-max`: maximum $\eta$
- `cone-radius`: radius used in Jet algorithm
- `obs`: observable[observable]
- `units`: `pb` for pico barn and `nb` for nano barn[unit]
- `value`: experimental values of observable
- `plot-factor`: multiplier for separating $\eta$ bins in plots

[observable]: `<` and `>` can only be used in pairs to represent averaging.
[unit]: Values of `units` have to be the same for the whole dataset, because the numeric unit conversion factor is read in only based on the first entry.

## uncertainties and corrections

When either the systematic or statistical uncertainties have a positive and negative and are different in magnitude, only the one with a larger magnitude is written in the data file.

- `_c` means correlated and and `_u` means uncorrelated

- `%` means the uncertainty or normalization or other types of corrections is a percentage

- `norm` is reserved for normalization

- `parton-to-hadron` is the correction from parton level calculation to hadron level, it is 1 plus the relative correction

## CDF
- We have examined carefully CJ15, the correlated systematic uncertainties are added in quadrature.
    - `./fitpack/theory/altpfit15.f`
        1. subroutine `fildat` instance in line 1857
        2. special treatment for CDF data in line 1947
    - `./fitpack/util/minim15.f`
        1. subroutine `DATSCN` in line 250
        2. special treatment for CDF data in line 312
- CJ15 has a covariance matrix for CDF run I data.
- [Joseph Owens](mailto:owens@hep.fsu.edu) claims that using the covariance matrix or adding everything in quadrature give same low $\chi ^2$.
- [Joseph Owens](mailto:owens@hep.fsu.edu) that for CDF run II data, CJ used the errors added in quadrature and the covariance matrix is for CDF run I data.
