/* Generated from /home/skimo/git/cloog-parma/cloog-core/test/vivien.cloog by CLooG 0.14.0-285-g341b1cd gmp bits. */
if (28*n >= -27) {
  for (p1=-54*n+4;p1<=4;p1++) {
    if (p1%2 == 0) {
      S1((p1-2)/2);
    }
  }
  if (n >= 1) {
    S3(1);
  }
  if ((n <= 1) && (n >= 2)) {
    S1(2);
  }
  for (p1=max(max(5,-54*n+4),4*n+2);p1<=6;p1++) {
    if (p1%2 == 0) {
      S1((p1-2)/2);
    }
  }
  if (n >= 2) {
    S4(1,2);
    S1(2);
    S6(1,2);
  }
  for (p1=7;p1<=min(9,4*n-2);p1++) {
    for (p2=ceild(-p1+2,4);p2<=min(-1,floord(-p1+2*n,2));p2++) {
      if (p1%2 == 0) {
        S4(-p2,(p1+2*p2)/2);
      }
    }
    if ((p1+3)%4 == 0) {
      S3((p1-1)/4);
    }
    if (p1%2 == 0) {
      S1((p1-2)/2);
    }
    if (p1 <= 2*n+2) {
      if (p1%2 == 0) {
        S6(1,(p1-2)/2);
      }
      if ((p1+1)%2 == 0) {
        S2((p1-3)/2,1);
      }
    }
    for (p2=max(ceild(-p1+2*n+5,2),ceild(p1-2*n-1,2));p2<=floord(p1-3,4);p2++) {
      if ((p1+1)%2 == 0) {
        S2((p1-2*p2-1)/2,p2);
      }
    }
  }
  for (p1=10;p1<=min(2*n+58,4*n-2);p1++) {
    for (p2=ceild(-p1+2,4);p2<=min(floord(-p1+2*n,2),floord(-p1+5,4));p2++) {
      if (p1%2 == 0) {
        S4(-p2,(p1+2*p2)/2);
      }
    }
    for (p2=ceild(-p1+6,4);p2<=min(min(-1,floord(-p1+2*n,2)),floord(-p1+9,4));p2++) {
      if (p1%2 == 0) {
        S4(-p2,(p1+2*p2)/2);
      }
      for (p3=1;p3<=-p2;p3++) {
        if (p1%2 == 0) {
          S5(-p2+1,(p1+2*p2-2)/2,p3);
        }
      }
    }
    for (p2=max(ceild(-p1+2*n+1,2),ceild(-p1+6,4));p2<=min(min(-1,floord(-p1+2*n+2,2)),floord(-p1+9,4));p2++) {
      for (p3=1;p3<=-p2;p3++) {
        if (p1%2 == 0) {
          S5(-p2+1,(p1+2*p2-2)/2,p3);
        }
      }
    }
    for (p2=ceild(-p1+10,4);p2<=min(-1,floord(-p1+2*n,2));p2++) {
      if (p1%2 == 0) {
        S4(-p2,(p1+2*p2)/2);
      }
      if (p1%2 == 0) {
        S6(-p2+2,(p1+2*p2-4)/2);
      }
      for (p3=1;p3<=-p2;p3++) {
        if (p1%2 == 0) {
          S5(-p2+1,(p1+2*p2-2)/2,p3);
        }
      }
    }
    for (p2=max(ceild(-p1+2*n+1,2),ceild(-p1+10,4));p2<=min(-1,floord(-p1+2*n+2,2));p2++) {
      if (p1%2 == 0) {
        S6(-p2+2,(p1+2*p2-4)/2);
      }
      for (p3=1;p3<=-p2;p3++) {
        if (p1%2 == 0) {
          S5(-p2+1,(p1+2*p2-2)/2,p3);
        }
      }
    }
    for (p2=max(ceild(-p1+2*n+3,2),ceild(-p1+10,4));p2<=min(-1,floord(-p1+2*n+4,2));p2++) {
      if (p1%2 == 0) {
        S6(-p2+2,(p1+2*p2-4)/2);
      }
    }
    if (p1 <= 2*n+4) {
      if (p1%2 == 0) {
        S6(2,(p1-4)/2);
      }
      if ((p1+3)%4 == 0) {
        S3((p1-1)/4);
      }
      if (p1%2 == 0) {
        S1((p1-2)/2);
      }
    }
    if (p1 >= 2*n+5) {
      if ((p1+3)%4 == 0) {
        S3((p1-1)/4);
      }
      if (p1%2 == 0) {
        S1((p1-2)/2);
      }
    }
    if (p1 <= 2*n+2) {
      if (p1%2 == 0) {
        S6(1,(p1-2)/2);
      }
      if ((p1+1)%2 == 0) {
        S2((p1-3)/2,1);
      }
    }
    for (p2=2;p2<=min(floord(-p1+2*n+4,2),floord(p1-3,4));p2++) {
      if ((p1+1)%2 == 0) {
        S2((p1-2*p2-1)/2,p2);
      }
    }
    for (p2=max(ceild(-p1+2*n+5,2),ceild(p1-2*n-1,2));p2<=floord(p1-3,4);p2++) {
      if ((p1+1)%2 == 0) {
        S2((p1-2*p2-1)/2,p2);
      }
    }
  }
  if ((n >= 2) && (n <= 29)) {
    S2(n,n-1);
  }
  for (p1=max(7,4*n);p1<=min(2*n+58,4*n+1);p1++) {
    if ((p1+3)%4 == 0) {
      S3((p1-1)/4);
    }
    if (p1%2 == 0) {
      S1((p1-2)/2);
    }
  }
  for (p1=max(max(7,-54*n+4),4*n+2);p1<=2*n+58;p1++) {
    if (p1%2 == 0) {
      S1((p1-2)/2);
    }
  }
  for (p1=2*n+59;p1<=4*n-2;p1++) {
    for (p2=ceild(-p1+2,4);p2<=min(floord(-p1+2*n,2),floord(-p1+5,4));p2++) {
      if (p1%2 == 0) {
        S4(-p2,(p1+2*p2)/2);
      }
    }
    for (p2=ceild(-p1+6,4);p2<=min(floord(-p1+2*n,2),floord(-p1+9,4));p2++) {
      if (p1%2 == 0) {
        S4(-p2,(p1+2*p2)/2);
      }
      for (p3=1;p3<=-p2;p3++) {
        if (p1%2 == 0) {
          S5(-p2+1,(p1+2*p2-2)/2,p3);
        }
      }
    }
    for (p2=max(ceild(-p1+2*n+1,2),ceild(-p1+6,4));p2<=min(floord(-p1+2*n+2,2),floord(-p1+9,4));p2++) {
      for (p3=1;p3<=-p2;p3++) {
        if (p1%2 == 0) {
          S5(-p2+1,(p1+2*p2-2)/2,p3);
        }
      }
    }
    for (p2=ceild(-p1+10,4);p2<=floord(-p1+2*n,2);p2++) {
      if (p1%2 == 0) {
        S4(-p2,(p1+2*p2)/2);
      }
      if (p1%2 == 0) {
        S6(-p2+2,(p1+2*p2-4)/2);
      }
      for (p3=1;p3<=-p2;p3++) {
        if (p1%2 == 0) {
          S5(-p2+1,(p1+2*p2-2)/2,p3);
        }
      }
    }
    for (p2=max(ceild(-p1+2*n+1,2),ceild(-p1+10,4));p2<=floord(-p1+2*n+2,2);p2++) {
      if (p1%2 == 0) {
        S6(-p2+2,(p1+2*p2-4)/2);
      }
      for (p3=1;p3<=-p2;p3++) {
        if (p1%2 == 0) {
          S5(-p2+1,(p1+2*p2-2)/2,p3);
        }
      }
    }
    for (p2=max(ceild(-p1+2*n+3,2),ceild(-p1+10,4));p2<=floord(-p1+2*n+4,2);p2++) {
      if (p1%2 == 0) {
        S6(-p2+2,(p1+2*p2-4)/2);
      }
    }
    if ((p1+3)%4 == 0) {
      S3((p1-1)/4);
    }
    for (p2=ceild(p1-2*n-1,2);p2<=floord(p1-3,4);p2++) {
      if ((p1+1)%2 == 0) {
        S2((p1-2*p2-1)/2,p2);
      }
    }
  }
  if (n >= 30) {
    S2(n,n-1);
  }
  for (p1=max(4*n,2*n+59);p1<=4*n+1;p1++) {
    if ((p1+3)%4 == 0) {
      S3((p1-1)/4);
    }
  }
}
