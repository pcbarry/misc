# unpolarized Jet data from Tevatron and RHIC

## data tables

| index | ref                    | normalization | collision | year                                   | collaboration    | status |
| ----- | -----                  | -----         | -----     | -----                                  | -----            | -----  |
| 20001 | [drupal][link.20001]   | no            | `pp`      | 2006 paper on 2003 and 2004 data       | STAR             | ready  |
| 20002 | [drupal][link.20002]   | yes           | `pp`      | 2012 paper on 2005 data                | STAR             | ready  |
| 20003 | [drupal][link.20003]   | yes           | `pp`      | 2012 paper on 2006 data                | STAR             | ready  |
| 20004 | [drupal][link.20004]   | yes           | `pp`      | 2015 paper on 2009 data                | STAR             | ready  |
| 20005 | [phenix][link.20005]   | no            | `pp`      | 2011 paper on 2005 data                | PHENIX           | ready  |
| 20006 | [drupal][link.20006]   | yes           | `pp`      | 2019 paper on 2012 data                | STAR             | ready  |

[link.20001]: https://drupal.star.bnl.gov/STAR/files/starpublications/68/data.html
[link.20002]: https://drupal.star.bnl.gov/STAR/files/starpublications/188/data.html
[link.20003]: https://drupal.star.bnl.gov/STAR/files/starpublications/188/data.html
[link.20004]: https://drupal.star.bnl.gov/STAR/files/starpublications/217/data.html
[link.20005]: https://www.phenix.bnl.gov/phenix/WWW/info/data/ppg093_data.html
[link.20006]: https://drupal.star.bnl.gov/STAR/files/starpublications/310/data.html

## observables

- $A_{L L}$

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

[observable]: `<` and `>` can only be used in pairs to represent averaging.
[unit]: Values of `units` have to be the same for the whole dataset, because the numeric unit conversion factor is read in only based on the first entry.

## uncertainties and corrections

When either the systematic or statistical uncertainties have a positive and negative and are different in magnitude, only the one with a larger magnitude is written in the data file.

- `_c` means correlated and and `_u` means uncorrelated

- `%` means the uncertainty or normalization or other types of corrections is a percentage

- `norm` is reserved for normalization

- `parton-to-hadron` is the correction from parton level calculation to hadron level, it is 1 plus the relative correction

## STAR 2019 paper on 2012 data
According to the paper, part of the systematic uncertainties in STAR 2019 paper on 2012 data are uncorrelated.
