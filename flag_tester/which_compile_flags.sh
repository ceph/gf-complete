if [ -n "$1" ]; then
  CC=$1
else
  CC=cc
fi

$CC flag_test.c -o flag_test 2> /dev/null
if [ -e "flag_test" ]; then
  OUTPUT=`./flag_test $CC 2> /dev/null`
  if [ -n "$OUTPUT" ]; then
    echo "$OUTPUT"
  else
    printf "CFLAGS = -O3\nLDFLAGS = -O3\n"
  fi
else
  printf "$CC failed to compile flag_test.c\n"
fi

rm sse4 sse2 ssse3 pclmul diff.txt flag_test temp.txt 2> /dev/null
