#!/usr/bin/env gawk -f

{
  i++;
  x[i]=$1;
  y[i]=$2;
}
END{
  n=i;
  s=0.0;
  for (i=1; i<=n-1; i++) {
    s=s+0.5*(y[i]+y[i+1])*(x[i+1]-x[i]);
    if (x[i]>18.5) { break; }
  }
  for (i=1; i<=n; i++) {
    print x[i],y[i]/s;
#    print x[i],y[i];
  }
}

