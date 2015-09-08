/* Generated from /home/skimo/git/cloog-parma/cloog-core/test/logopar.cloog by CLooG 0.14.0-285-g341b1cd gmp bits. */
for (j=0;j<=m;j++) {
  S1(1,j);
}
for (i=2;i<=n;i++) {
  for (j=0;j<=i-2;j++) {
    S2(i,j);
  }
  for (j=i-1;j<=n;j++) {
    S1(i,j);
    S2(i,j);
  }
  for (j=n+1;j<=m;j++) {
    S1(i,j);
  }
}
for (i=n+1;i<=m+1;i++) {
  for (j=i-1;j<=m;j++) {
    S1(i,j);
  }
}
