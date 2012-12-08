($1 == "Seed:") { l = 0; n++; t=0 }
{ if (l >= 1 && l <= 4) {
    t += $4
    if (l == 4) avg = t/4.0
  }
  if (l == 5) {
    printf("xaxis max %d hash_label at %d : %s\n", n+1, n, $0 )
    printf("newcurve marktype xbar marksize 1 cfill 1 1 0 pts %d %.2lf\n", n, avg);
  }
  l++
}
