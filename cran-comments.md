# Cran comments

This is a binary resubmit for `roxde2parse` (since dparser was just accepted).  Since it generates the
code for 'nlmixr2'/'rxode2' in C this is also tested against the
development versions to make sure that:

-Wstrict-prototypes

There were a few minor errors here too, which were fixed.

This version also adds support for the upcoming `nlmixr2random` which
splits `rxode2` further to hopefully have less than a 10 minute
install and check for the `rxode2` package.

