/* Generated from /home/skimo/git/cloog-parma/cloog-core/test/darte.cloog by CLooG 0.14.0-285-g341b1cd gmp bits. */
if (n >= 1) {
  for (t3=n+3;t3<=3*n+1;t3++) {
    if ((t3+n+1)%2 == 0) {
      S1(1,n,(t3-n-1)/2);
    }
  }
  for (t1=-n+2;t1<=n-1;t1++) {
    for (t2=max(-t1+2,t1+2);t2<=min(-t1+4,-t1+2*n);t2++) {
      for (t3=t2+2;t3<=t2+2*n;t3++) {
        if ((t1+t2)%2 == 0) {
          if ((t1+t3)%2 == 0) {
            S1((t1+t2)/2,(-t1+t2)/2,(-t2+t3)/2);
          }
        }
      }
    }
    if (t1 >= 2) {
      for (t3=t1+4;t3<=t1+2*n+2;t3++) {
        if ((t1+t3)%2 == 0) {
          S1(t1+1,1,(-t1+t3-2)/2);
        }
      }
    }
    for (t2=max(-t1+5,t1+3);t2<=min(-t1+2*n,t1+2*n);t2++) {
      for (t3=1;t3<=min(n,t2+1);t3++) {
        if ((t1+t2+1)%2 == 0) {
          S2((t1+t2-3)/2,(-t1+t2-1)/2,t3);
        }
      }
      for (t3=t2+2;t3<=n;t3++) {
        if ((t1+t2+1)%2 == 0) {
          S2((t1+t2-3)/2,(-t1+t2-1)/2,t3);
        }
        if ((t1+t2)%2 == 0) {
          if ((t1+t3)%2 == 0) {
            S1((t1+t2)/2,(-t1+t2)/2,(-t2+t3)/2);
          }
        }
      }
      for (t3=max(n+1,t2+2);t3<=t2+2*n;t3++) {
        if ((t1+t2)%2 == 0) {
          if ((t1+t3)%2 == 0) {
            S1((t1+t2)/2,(-t1+t2)/2,(-t2+t3)/2);
          }
        }
      }
    }
    for (t2=max(-t1+5,-t1+2*n+1);t2<=min(-t1+2*n+3,t1+2*n+1);t2++) {
      for (t3=1;t3<=n;t3++) {
        if ((t1+t2+1)%2 == 0) {
          S2((t1+t2-3)/2,(-t1+t2-1)/2,t3);
        }
      }
    }
    if (t1 <= -1) {
      for (t3=1;t3<=n;t3++) {
        S2(t1+n-1,n,t3);
      }
    }
  }
  for (t3=1;t3<=n;t3++) {
    S2(n,1,t3);
  }
}
