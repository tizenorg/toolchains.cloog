/* Generated from ../../../git/cloog/test/tiling.cloog by CLooG 0.14.0-253-ge300ff5 gmp bits in 0.00s. */
for (ii=0;ii<=floord(n,10);ii++) {
  for (i=max(0,10*ii);i<=min(n,10*ii+9);i++) {
    S1(ii,i);
  }
}