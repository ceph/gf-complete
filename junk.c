#include <stdio.h>

main()
{
  int size, iterations;
  double ds, di, elapsed;

  elapsed = 0.614553;
  size = 8192;
  iterations = 655360;

  ds = size;
  di = iterations;

  printf("%10.3lf\n", ((double) (size*iterations)) / (1024 * 1024 * elapsed));
  printf("%10.3lf\n", ds * di / 1024.0 / 1024.0 / elapsed);
}
