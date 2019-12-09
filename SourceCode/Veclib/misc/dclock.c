/*
 * Double-precision timing routine
 */

#include <time.h>

#if !defined(i860) && !defined(dclock)

double dclock(void)
{
  static double tps = 1.0 / CLOCKS_PER_SEC;
  return (double) clock() * tps;
}

#endif
