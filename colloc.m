function [r, A, B, q]=colloc(n,left,right)
% colloc: Calculate collocation weights
%         [r, A, B, q] = colloc( n [,'left'] [,'right'])
% inputs:
%               n - number of interior node points
%          'left' - include left boundary
%         'right' - include right bounary also
% outputs:
%          r - vector of roots
%          A - Matrix of first derivative weights
%          B - Matrix of second derivative weights
%          q - Quadrature weights.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 1996, 1997 John W. Eaton
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, write to the Free
% Software Foundation, 59 Temple Place - Suite 330, Boston, MA
% 02111-1307, USA.
%
% Adapted from Octave's colloc.cc by Steve Swinnea.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n0 = 0 ; n1 = 0;
if (nargin > 1)
  if (strcmp(left,'left') | strcmp(left,'l') )
    n0 = 1;
  elseif (left == 0 | left == 1 )
    n0 = left;
  else
    error('Second argument should be the string left or l')
  end
end
if (nargin > 2)
  if (strcmp(right,'right') | strcmp(right,'r') )
    n1 = 1;
  elseif ( right == 1 | right == 0 )
    n1 = right;  
  else
    error('Third argument should be the string right or r')
  end
end

[dif1,dif2,dif3,r]=jcobi(n,n0,n1,0,0);
q = dfopr(n,n0,n1,0,3,dif1,dif2,dif3,r);
for i=1:(n+n0+n1)
    vect = dfopr(n,n0,n1,i,1,dif1,dif2,dif3,r);
    A(i,:) = vect';
end
for i=1:(n+n0+n1)
    vect = dfopr(n,n0,n1,i,2,dif1,dif2,dif3,r);
    B(i,:) = vect';
end

%%%%%%  jcobi  %%%%%%%
function [dif1,dif2,dif3,root]=jcobi(n,n0,n1,alpha,beta)
if (n0 ~= 0) & (n0 ~= 1)
    error('** VILERR : Illegal value %  N0 ');
end
if (n1 ~= 0) & (n1 ~= 1)
    error('** VILERR : Illegal value for N1 ');
end
if (n+n0+n1 < 1)
    error('** VILERR : Number of interpolation points less than 1');
end
%
% -- FIRST EVALUATION OF COEFFICIENTS IN RECURSION FORMULAS.
% -- RECURSION COEFFICIENTS ARE STORED IN DIF1 AND DIF2.
%
nt = n+n0+n1;
dif1=zeros(nt,1);
dif2=zeros(nt,1);
dif3=zeros(nt,1);
root=zeros(nt,1);
ab = alpha+beta;
ad = beta-alpha;
ap = beta*alpha;
dif1(1) = (ad/(ab+2)+1)/2;
dif2(1) = 0;

if (n >= 2)
    for i=2:n
        z1 = i-1;
        z = ab + 2*z1;
        dif1(i) = (ab*ad/z/(z+2)+1)/2;
        if (i == 2 )
            dif2(i) = (ab+ap+z1)/z/z/(z+1);
        else
            z = z*z;
            y = z1*(ab+z1);
            y = y*(ap+y);
            dif2(i) = y/z/(z-1);
        end
    end
end
%
% -- ROOT DETERMINATION BY NEWTON METHOD WITH SUPPRESSION OF
% -- PREVIOUSLY DETERMINED ROOTS
%
x = 0;
for i=1:n
  z = 1;
  while ( abs(z) > 1e-9 )
    xd = 0;
    xn = 1;
    xd1 = 0;
    xn1 = 0;
    for j=1:n
        xp = (dif1(j)-x)*xn - dif2(j)*xd;
        xp1 = (dif1(j)-x)*xn1 - dif2(j)*xd1 - xn;
        xd = xn;
        xd1 = xn1;
        xn = xp;
        xn1 = xp1;
    end
    zc = 1;
    z = xn/xn1;
    if ( i ~= 1 )
        for j = 2:i
            zc = zc - z/(x-root(j-1));
        end
    end
    z = z/zc;
    x = x-z;
  end
  root(i) = x;
  x = x +.0001;
end
% 
% -- ADD INTERPOLATION POINTS AT X = 0 AND/OR X = 1
% 
if (n0 ~= 0)
    root = [ 0 ; root(1:nt-1) ];
end
if (n1 == 1)
    root(nt) = 1;
end
[dif1 dif2 dif3] = dif( root );

%%%%% dfopr %%%%%%
function vect = dfopr(n,n0,n1,i,id,dif1,dif2,dif3,root)
nt = n+n0+n1;
vect = zeros(nt,1);
if (n0 ~= 0) & (n0 ~= 1)
    error('** VILERR : Illegal value %  N0 ');
end
if (n1 ~= 0) & (n1 ~= 1)
    error('** VILERR : Illegal value for N1 ');
end
if (nt < 1)
    error('** VILERR : Number of interpolation points less than 1');
end
if (id ~= 1 & id ~= 2 & id ~= 3 )
    error('** VILERR : Illegal ID in DFOPR ')
end
if ( id ~= 3 )
    if ( i < 1 )
        error('** VILERR : Index less than zero in DFOPR ')
    end
    if ( i > nt )
        error('** VILERR : Index greater than NTOTAL in DFOPR ')
    end
end

% 
% -- EVALUATE DISCRETIZATION MATRICES AND GAUSSIAN QUADRATURE
% -- WEIGHTS.  QUADRATURE WEIGHTS ARE NORMALIZED TO SUM TO ONE.
% 
if ( id ~= 3 )
    for j = 1:nt
        if (j == i)
            if (id == 1)
                vect(i) = dif2(i)/dif1(i)/2;
            else
                vect(i) = dif3(i)/dif1(i)/3;
            end
        else
            y = root(i)-root(j);
            vect(j) = dif1(i)/dif1(j)/y;
            if (id == 2 )
                vect(j)=vect(j)*(dif2(i)/dif1(i)-2/y);
            end
        end
    end
else
    y=0;
    for j = 1:nt
        x = root(j);
        ax = x*(1-x);
        if (n0 == 0)
            ax = ax/x/x;
        end
        if (n1 == 0)
            ax = ax/(1-x)/(1-x);
        end
        vect(j) = ax/dif1(j)^2;
        y = y + vect(j);
    end
    vect = vect/y;
end

%%%%% dif %%%%%
function [dif1,dif2,dif3] = dif( root )
nt = length( root );
dif1 = zeros(nt,1);
dif2 = zeros(nt,1);
dif3 = zeros(nt,1);
if ( nt < 1 )
    error('** VILERR : Number of interpolation points less than 1');
end
for i = 1:nt
    x = root(i);
    dif1(i) = 1;
    dif2(i) = 0;
    dif3(i) = 0;
    for j = 1:nt
        if ( j ~= i)
            y = x - root(j);
            dif3(i) = y*dif3(i) + 3*dif2(i);
            dif2(i) = y*dif2(i) + 2*dif1(i);
            dif1(i) = y*dif1(i);
        end
    end
end

