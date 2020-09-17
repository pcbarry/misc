# pion-tungsten (per nucleon) Drell-Yan cross sections from E615 experiment at FNAL.
  
## Data sets:
* http://hepdata.cedar.ac.uk/review/dy/e615.shtml

## JAM database (converted into xlsx from ASCII files in CJ database)

The table below indexes and summarizes the E610 pion-tungsten Drell-Yan cross sections.

### Observables

* dsigma/dsqrt(tau)dx
* dsigma/dpTdQ
* dsigma/dpTdxF

### Columns:

##File 1000
- tau_min  = tau_min [GeV]
- tau_max  = tau_max [GeV]
- x_min  = xF_min [GeV]
- x_max  = xF_max [GeV]
- dsigma/dsqrt(tau)dx = dsigma/dsqrt(tau)dxF 
- stat_error  = statistical uncertainty

##File 1001
- m_min  = Q_min [GeV]
- m_max  = Q_max [GeV]
- pT     = pT
- dsigma/dpTdm = dsigma/dpTdQ
- stat_error  = statistical uncertainty

##File 1002
- x_min  = xF_min [GeV]
- x_max  = xF_max [GeV]
- pT     = pT
- dsigma/dpTdx = dsigma/dpTdxF
- stat_error  = statistical uncertainty

## Data table

| index | ref              | process | target | obs                 | experiment    | status |
| :--:  | :--:             | :--:    | :--:   | :--:                | :--:          | :--:   |
| 1000  | [link][ref1000]  | DY      |tungsten| dsigma/dsqrt(tau)dx | Fermilab E615 | FINAL  |
| 1001  | [link][ref1001]  | DY      |tungsten| dsigma/dpTdQ        | Fermilab E615 | FINAL  |
| 1002  | [link][ref1002]  | DY      |tungsten| dsigma/dpTdxF       | Fermilab E615 | FINAL  |

[ref1000]: https://journals.aps.org/prd/pdf/10.1103/PhysRevD.39.92
[ref1001]: https://journals.aps.org/prd/pdf/10.1103/PhysRevD.39.92
[ref1002]: https://journals.aps.org/prd/pdf/10.1103/PhysRevD.39.92
