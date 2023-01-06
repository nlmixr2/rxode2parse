
/*
strncmpci.c

- A 'c'ase 'i'nsensitive version of `strncmp()`.
- See references below for more info., including documentation for `strncmp()`, as well as my
Stack Overflow answer where I present my `strncmpci()` function below.

Gabriel Staples
www.ElectricRCAircraftGuy.com
Written: 21 Mar. 2019
Updated: 21 Oct. 2020
- moved to this git repo; see `git log` history after that

Matthew Fidler, Incorporated in rxode2parse AND removed test code here (will adapt to a R test case)

- A R interface for testing


References:
1. [my own answer (Gabriel Staples)] https://stackoverflow.com/questions/5820810/case-insensitive-string-comp-in-c/55293507#55293507
2. https://en.cppreference.com/w/cpp/string/byte/strncmp
3. http://www.cplusplus.com/reference/cstring/strncmp/

STATUS:
IT WORKS! ALL UNIT TESTS PASS!

*/

#define USE_FC_LEN_T
#define STRICT_R_HEADERS
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#include <assert.h>
#include <stdbool.h>
#include <ctype.h> // for `tolower()`
#include <limits.h> // for `INT_MIN`
#include <stdio.h>
#include <string.h>
#include "strncmpi.h"

typedef struct data_s
{
    int error_count;
} data_t;

// Data struct used to safely contain and pass around global data
data_t globals = {
    .error_count = 0,
};

// TODO: Make a version of this code which also works on Unicode's UTF-8 implementation (character
// encoding)! Add it to my answer here too: https://stackoverflow.com/a/55293507/4561887.

/// \brief      Perform a case-insensitive string compare (`strncmp()` case-insensitive) to see
///             if two C-strings are equal.
/// \note       1. Identical to `strncmp()` except:
///               1. It is case-insensitive.
///               2. The behavior is NOT undefined (it is well-defined) if either string is a null
///               ptr. Regular `strncmp()` has undefined behavior if either string is a null ptr
///               (see: https://en.cppreference.com/w/cpp/string/byte/strncmp).
///               3. It returns `INT_MIN` as a special sentinel value for certain errors.
///             - Posted as an answer here: https://stackoverflow.com/a/55293507/4561887.
///               - Aided/inspired, in part, by `strcicmp()` here:
///                 https://stackoverflow.com/a/5820991/4561887.
/// \param[in]  str1        C string 1 to be compared.
/// \param[in]  str2        C string 2 to be compared.
/// \param[in]  num         max number of chars to compare
/// \return     A comparison code (identical to `strncmp()`, except with the addition
///             of `INT_MIN` as a special sentinel value):
///
///             INT_MIN (usually -2147483648 for int32_t integers)  Invalid arguments (one or both
///                      of the input strings is a NULL pointer).
///             <0       The first character that does not match has a lower value in str1 than
///                      in str2.
///              0       The contents of both strings are equal.
///             >0       The first character that does not match has a greater value in str1 than
///                      in str2.
int strncmpci(const char * str1, const char * str2, size_t num)
{
    int ret_code = 0;
    size_t chars_compared = 0;

    // Check for NULL pointers
    if (!str1 || !str2)
    {
        ret_code = INT_MIN;
        return ret_code;
    }

    // Continue doing case-insensitive comparisons, one-character-at-a-time, of `str1` to `str2`, so
    // long as 1st: we have not yet compared the requested number of chars, and 2nd: the next char
    // of at least *one* of the strings is not zero (the null terminator for a C-string), meaning
    // that string still has more characters in it.
    // Note: you MUST check `(chars_compared < num)` FIRST or else dereferencing (reading) `str1` or
    // `str2` via `*str1` and `*str2`, respectively, is undefined behavior if you are reading one or
    // both of these C-strings outside of their array bounds.
    while ((chars_compared < num) && (*str1 || *str2))
    {
        ret_code = tolower((int)(*str1)) - tolower((int)(*str2));
        if (ret_code != 0)
        {
            // The 2 chars just compared don't match
            break;
        }
        chars_compared++;
        str1++;
        str2++;
    }

    return ret_code;
}

extern int rxstrcmpi(const char * str1, const char * str2) {
  return strncmpci(str1, str2, INT_MAX);  
}


SEXP _rxode2_parse_strncmpci(SEXP str1, SEXP str2, SEXP num) {
  SEXP reti = PROTECT(Rf_allocVector(INTSXP, 1));
  INTEGER(reti)[0] = strncmpci(R_CHAR(STRING_ELT(str1, 0)),
                               R_CHAR(STRING_ELT(str2, 0)),
                               INTEGER(num)[0]);
  UNPROTECT(1);
  return reti;
}
