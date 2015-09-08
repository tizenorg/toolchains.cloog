/* Generated from /home/skimo/git/cloog-parma/cloog-core/test/thomasset.cloog by CLooG 0.14.0-285-g341b1cd gmp bits. */
if (n >= 1) {
  for (c1=0;c1<=floord(n-5,3);c1++) {
    for (i=max(1,3*c1+1);i<=3*c1+3;i++) {
      S1(i,c1);
    }
  }
  for (c1=max(0,ceild(n-4,3));c1<=floord(n-1,3);c1++) {
    if (c1 <= 0) {
      S1(1,c1);
      for (j=1;j<=min(n,3*c1-n+5);j++) {
        for (k=0;k<=min(0,floord(3*c1-j-n+4,3));k++) {
          for (p=max(ceild(n-2,3),ceild(3*c1-j-3*k,3));p<=min(floord(n,3),floord(3*c1-j-3*k+2,3));p++) {
            S2(1,j,k,p,c1-k-p);
          }
        }
      }
    }
    if (c1 >= 1) {
      for (j=1;j<=min(n,3*c1-n+5);j++) {
        for (k=0;k<=min(0,floord(3*c1-j-n+4,3));k++) {
          for (p=max(ceild(n-2,3),ceild(3*c1-j-3*k,3));p<=min(floord(n,3),floord(3*c1-j-3*k+2,3));p++) {
            S2(1,j,k,p,c1-k-p);
          }
        }
      }
    }
    for (i=max(2,3*c1+1);i<=min(n,3*c1+3);i++) {
      S1(i,c1);
    }
    for (c2=1;c2<=n-1;c2++) {
      for (j=1;j<=min(n,3*c1-n+5);j++) {
        for (k=0;k<=min(0,floord(3*c1-j-n+4,3));k++) {
          for (p=max(ceild(n-2,3),ceild(3*c1-j-3*k,3));p<=min(floord(n,3),floord(3*c1-j-3*k+2,3));p++) {
            S2(c2+1,j,k,p,c1-k-p);
          }
        }
      }
    }
  }
  for (c1=ceild(n,3);c1<=floord(2*n+1,3);c1++) {
    for (c2=0;c2<=n-1;c2++) {
      for (j=max(1,3*c1-n-1);j<=min(n,3*c1-n+5);j++) {
        for (k=max(0,ceild(3*c1-j-n,3));k<=min(0,floord(3*c1-j-n+4,3));k++) {
          for (p=max(ceild(n-2,3),ceild(3*c1-j-3*k,3));p<=min(floord(n,3),floord(3*c1-j-3*k+2,3));p++) {
            S2(c2+1,j,k,p,c1-k-p);
          }
        }
      }
    }
  }
}
