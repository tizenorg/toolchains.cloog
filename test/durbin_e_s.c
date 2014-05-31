/* Generated from /home/skimo/git/cloog-parma/cloog-core/test/durbin_e_s.cloog by CLooG 0.14.0-285-g341b1cd gmp bits. */
S4(1,0,0);
S7(1,0,0);
S8(1,0,3);
for (i=2;i<=9;i++) {
  if (i == 2) {
    S2(i,-7,0);
    S3(i,-7,1);
    S6(i,-7,2);
  }
  if (i >= 3) {
    S2(i,-7,0);
    S3(i,-7,1);
  }
  for (j=-6;j<=i-10;j++) {
    S3(i,j,1);
  }
  if (i == 9) {
    S3(i,0,1);
    S6(i,0,2);
    S8(i,0,3);
  }
  if ((i >= 3) && (i <= 8)) {
    S3(i,i-9,1);
    S6(i,i-9,2);
  }
  if (i <= 8) {
    S8(i,0,3);
  }
  for (j=1;j<=i-1;j++) {
    S5(i,j,3);
  }
}
S2(10,-7,0);
S3(10,-7,1);
for (j=-6;j<=0;j++) {
  S3(10,j,1);
}
S3(10,1,1);
S6(10,1,2);
S5(10,1,3);
S1(10,1,4);
for (j=2;j<=9;j++) {
  S5(10,j,3);
  S1(10,j,4);
}
S1(10,10,4);
