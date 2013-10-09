Run which_compile_flags.sh and it will print out the compile flags to use in
  GNUmakefile. By default, this script uses "cc" as its compiler but you can
  pass in the name of your compiler as an argument.

EXAMPLE: "./which_compile_flags.sh clang"

This script will run "clang" in the above example so be warned that if you type
something like "rm" for that argument, you get what you asked for.  Also, make
sure that the compiler that you pass to which_compile_flags.sh is the same as
the compiler in GNUmakefile.
