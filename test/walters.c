/* Generated from ../../../git/cloog/test/walters.cloog by CLooG 0.14.0-261-g26b30f1 gmp bits in 0.02s. */
div36 = 0;
S2(1,div36,1,0);
S4(1,div36,1,0);
S3(2,0,1,1);
S4(2,0,1,1);
for (i=3;i<=10;i++) {
  if ((i+2)%3 <= 1) {
    div36 = floord(i-1,3);
    if ((i+1)%3 <= 1) {
      div37 = floord(i+1,3);
      if (i%3 <= 1) {
        div38 = floord(i,3);
        S4(i,div36,div37,div38);
      }
      if ((i+1)%3 == 0) {
        S3(i,div36,div37,(i+1)/3);
        S4(i,div36,div37,(i+1)/3);
      }
    }
    if ((i+2)%3 == 0) {
      div38 = floord(i,3);
      S2(i,div36,(i+2)/3,div38);
      S4(i,div36,(i+2)/3,div38);
    }
  }
  if (i%3 == 0) {
    div37 = floord(i+1,3);
    div38 = floord(i,3);
    S1(i,i/3,div37,div38);
    S4(i,i/3,div37,div38);
  }
}
