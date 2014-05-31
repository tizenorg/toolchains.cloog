/* Generated from /home/skimo/git/cloog-parma/cloog-core/test/mode.cloog by CLooG 0.14.0-285-g341b1cd gmp bits. */
if (M >= 0) {
  if (N >= 0) {
    for (i=0;i<=M;i++) {
      for (j=0;j<=min(N,i);j++) {
        S1(i,j);
        S2(i,j);
      }
      for (j=N+1;j<=i;j++) {
        S1(i,j);
      }
      for (j=i+1;j<=N;j++) {
        S2(i,j);
      }
    }
  }
  if (N <= -1) {
    for (i=0;i<=M;i++) {
      for (j=0;j<=i;j++) {
        S1(i,j);
      }
    }
  }
}
