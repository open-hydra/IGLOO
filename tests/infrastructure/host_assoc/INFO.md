# INFO — host_assoc (internal procedures address the main program's arrays, unit)

`test_host_assoc` (ctest `test_host_assoc`, labels `unit;infrastructure;build;t1`).

## What it pins
A test program's internal procedures read and write the main program's variables by host association. With the
Intel compilers' OpenMP option, local variables are automatic: a main-program array of explicit shape that the main
program hands to an assumed-shape dummy (a descriptor) then lives on the main program's stack, while an internal
procedure addresses a static copy of it. The internal procedure reads zeros, or the type's default initialization,
and what it writes never reaches the main program. Scalars, allocatable arrays and arrays passed to explicit-shape
dummies are not affected, and neither are the locals of module procedures.

`tests/CMakeLists.txt` compiles every test executable with `-auto-scalar` under the Intel compilers (the compilers'
storage model without OpenMP: scalars automatic, arrays and derived types static; local allocatables are still
deallocated at return). This program fails when a test executable is built without it. gfortran keeps one copy of
the array and passes with or without the option.

## Tests
The main program fills `tab(1:10) = 1..10` and `recs(1:10)%a = 1..10` (`rec_t` initializes `a` to −1) and passes
both to assumed-shape dummies of a module procedure that leaves them unchanged.

| id | checks |
|---|---|
| H1 | an internal procedure reads `tab(i) = i` for all ten elements |
| H2 | after an internal procedure writes `tab(i) = −i`, the main program reads −i for all ten |
| H3 | an internal procedure reads `recs(i)%a = i` for all ten |

## RED (the option removed)

| build | H1 | H2 | H3 |
|---|---|---|---|
| ifx 2023.1 `-fiopenmp` | 10 of 10 differ (`tab(5)` reads 0) | 10 of 10 differ (`tab(5)` reads 5) | 10 of 10 differ (`recs(5)%a` reads −1) |
| ifx 2023.1 `-fiopenmp -auto-scalar` | pass | pass | pass |

The result does not depend on inlining: the same split shows at `-O0` and with `-auto` alone.
