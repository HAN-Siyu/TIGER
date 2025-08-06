# TIGERr 1.0.2 (06-Aug-2025)

* Added `set_seed` option to the `run_TIGER()` function for reproducible results.
* Pearson correlation-based feature selection now supports parallel computing to speed up large datasets and prevent memory issues caused by `cor()`.
* Now notifies when a metabolite lacks enough observations to compute correlation during feature selection.

# TIGERr 1.0.1 (01-Nov-2023)

* Support passing values to the argument `use` in  `cor()` to control how to deal with missing values when computing correlation.

# TIGERr 1.0.0 (04-Jan-2022)

* Updated reference information.
* Fixed the issues caused by dependencies.

# TIGERr 0.1.0 (02-Sep-2021)

* Original release
