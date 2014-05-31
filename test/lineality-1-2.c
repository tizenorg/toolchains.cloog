/* Generated from /home/skimo/git/cloog-parma/cloog-core/test/lineality-1-2.cloog by CLooG 0.14.0-285-g341b1cd gmp bits. */
for (i=1;i<=M;i++) {
  for (j=1;j<=i-1;j++) {
    S1(i,j);
  }
  S1(i,i);
  S2(i,i);
  for (j=i+1;j<=M;j++) {
    S1(i,j);
  }
}
