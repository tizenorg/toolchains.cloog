/* Generated from ../../../git/cloog/test/walters2.cloog by CLooG 0.14.0-249-g82c3b0f gmp bits in 0.01s. */
for (i=0;i<=51;i++) {
  S2(0,i);
}
for (j=1;j<=24;j++) {
  S2(j,0);
  for (i=1;i<=50;i++) {
    S1(j,i);
  }
  S2(j,51);
}
S2(25,0);
for (i=1;i<=51;i++) {
  S2(25,i);
}
