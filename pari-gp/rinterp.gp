/*
  Adapted from: http://axiom-wiki.newsynthesis.org/RationalInterpolation
*/
rinterp(xlist, ylist, m, k)=
{
  my(collist, resspace);
  if(length(xlist)!=length(ylist),
    error("Different number of points and values."));
  if(length(xlist)!=m+k+1,
    error("Wrong number of points."));
    
  collist = matrix(m+k+1, m+k+2, i, j, 1);
  for(j=1, max(m,k),
    for(i=1, m+k+1, collist[i,j+1] = collist[i,j]*xlist[i]));
  for(j=1, k+1,
    for(i=1, m+k+1, collist[i, m+k+3-j] = -collist[i,k+2-j]*ylist[i]));
  print(collist);
  resspace = matker(collist);
  print (resspace);
  Pol(vector(m+1,j,resspace[m+2-j,1]))/Pol(vector(k+1,j,resspace[m+k+3-j,1]))
};
addhelp(rinterp,"rinterp(x,y,m,k): rational interpolation for a set of m+k+1 points exactly. Both vectors x and y must have a length equal to m+k+1; x and y contain the coordinates of the m+k+1 points. Degree of the numerator will be m; degree of the denominator will be k.");
