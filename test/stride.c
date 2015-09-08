/* Generated from ../../..//git/cloog-parma/cloog-core/test/stride.cloog by CLooG 0.14.0-275-g15c65cd gmp bits. */
for (c1=3;c1<=100;c1++) {
  for (c2=max(ceild(4*c1-100,9),ceild(-c1+25,22));c2<=floord(c1,3);c2++) {
    if ((c1 == 25) && (c2 == 0)) {
      S1(25);
    }
    if (c1 == 3*c2) {
      if (c1%3 == 0) {
        S2(c1,c1/3);
      }
    }
  }
}
