# rxode2parse 2.0.14

* 'linCmt()' translations of 'alpha', 'beta', 'gamma', 'k21', 'k31',
  'vc' now error instead of ignoring 'gamma' and 'k31' to give 2 cmt
  solution

* transit compartment internal code now changes dose to 0.0 when no
  dose has been administered to the depot compartment. This way dosing
  to the central compartment (without dosing to the transit
  compartment) will not give a `NA` for the depot compartment (and
  consequently for the central compartment)

* Moved `rxDerived` here and added tests for it here as well.

* Moved `etTransParse` here and added tests for it here as well (makes
  up most of `etTrans`). In addition the following changes were made
  to `etTransParse()`/`etTrans()`:

  * The internal translation (`etTrans()`) will not drop times when
    infusions stop. Before, if the infusion stopped after the last
    observation the time when the infusion stopped would be dropped.
    This interferes with `linCmt()` models.

  * Breaking change/bug fix `evid=2` are considered observations when
    translating data to internal `rxode2` event structure

  * Fix edge case to find infusion duration when it is the first item
    of the dosing record at time 0.

 * Fixed a bug for certain infusions where the `rate`, `ii` and/or
   `ss` data items were dropped from the output when `addDosing=TRUE`


* Also have internal functions to convert between classic NONMEM
  events and rxode2 events

* Have an internal function that gives information on the linear
  compartmental model translation type, which could be useful for
  babelmixr2

* 'time' in model is now case insensitive

* Use function declaration in `rxode2parseGetTranslation()` to
  determine thread safety of functions available to rxode2

* Add check for correct number of function arguments to parser.

* Like R, known functions can be assigned as a variable and the
  function can still be called (while not changing the variable
  value).  For example you can have a variable `gamma` as well as a
  function `gamma()`.

* Fix garbled error messages that occur with certain messages.

* Fixed errors that occurred when using capitalized AMT variables in
  the model.

# rxode2parse 2.0.13

* Version bump for dparser (so binaries will be built correctly)

# rxode2parse 2.0.12

* Bug fix for strict prototypes

* Removed `sprintf` as noted by CRAN

* Made `rxode2parse` dll binary independent of `rxode2()`

# rxode2parse 2.0.11

* Bug fix for strict aliasing as requested by CRAN

# rxode2parse 2.0.10

* Use strict aliasing as requested by CRAN

# rxode2parse 2.0.9

* Initial release to split of rxode2parse from rxode2 to reduce
  compilation time of 'rxode2'

* Added a `NEWS.md` file to track changes to the package.
