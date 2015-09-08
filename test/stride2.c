/* Generated from ../../../git/cloog-parma/cloog-core/test/stride2.cloog by CLooG 0.14.0-275-g15c65cd gmp bits. */
for (c1=3;c1<=100;c1++) {
  for (c2=max(ceild(-c1+27,24),ceild(100*c1-2700,219));c2<=floord(c1,3);c2++) {
    if ((c1 == 27) && (c2 == 0)) {
      S1(27);
    }
    if (c1 == 3*c2) {
      if (c1%3 == 0) {
        S2(c1,c1/3);
      }
    }
  }
}
