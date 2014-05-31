/* Generated from /home/skimo/git/cloog-parma/cloog-core/test/./non_optimal/youcef.cloog by CLooG 0.14.0-285-g341b1cd gmp bits. */
for (i=0;i<=5;i++) {
  if (i == 5) {
    S1(i,5);
    S2(i,5);
    S3(i,5);
  }
  if (i <= 4) {
    S1(i,i);
    S2(i,i);
  }
  for (j=i+1;j<=4;j++) {
    S2(i,j);
  }
  if (i <= 4) {
    S2(i,5);
    S3(i,5);
  }
}
