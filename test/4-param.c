/* Generated from ../../../git/cloog/test/4-param.cloog by CLooG 0.14.0-249-g82c3b0f gmp bits in 0.01s. */
for (i=m;i<=min(min(n,q),p-1);i++) {
  S1(i);
}
for (i=p;i<=min(min(n,q),m-1);i++) {
  S2(i);
}
for (i=max(m,p);i<=min(n,q);i++) {
  S1(i);
  S2(i);
}
for (i=max(m,q+1);i<=n;i++) {
  S1(i);
}
for (i=max(p,n+1);i<=q;i++) {
  S2(i);
}
