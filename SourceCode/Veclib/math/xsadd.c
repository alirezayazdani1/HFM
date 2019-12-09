/*
 * scalar plus a vector
 */

void dsadd (int n, double alpha, double *x, int incx, double *y, int incy)
{
  while (n--) {
    *y = alpha + *x;
    x += incx;
    y += incy;
  }
  return;
}
