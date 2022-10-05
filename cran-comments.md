# Cran comments

While this was *just* accepted, Prof Riply sent the following email:

While it applied to 'rxode2', the portion that applies to 'rxode2' has moved here.

Here is the email below:

[Packages listed above] import headers from
the one on the left with some 'prototypes' fn () .  Compiling using

-Wstrict-prototypes

in CFLAGS (with gcc or clang < 15: this is implied by -pedantic in clang
15) shows warnings that are starting to show up on the CRAN results
pages, and will be used for CRAN incoming checks.

Some of these seem intended to be a function with no arguments, for
which the prototype is fn (void) .  But for others, consistency checks
are being circumvented by not using a full prototype.  In packages not
reported here these have shown up several bugs already ....

Please correct (in the packages before the colon) before 2022-11-04.

* This fix corrects this for rxode2parse
