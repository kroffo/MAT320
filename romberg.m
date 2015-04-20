function [ ans, flag, iterations ] = romberg(f, a, b, tol, max)
% f is the function to be integrated on the interval a to b as defined
% by the input to the function. tol is the tolerance within which
% successive romberg diagonal differences must be in order for the method
% to be successful. max is the maximum number of iterations, or the maximum
% number of diagonals to calculate.
%
% The function returns ans, the approximated value of the definite integral,
% flag, which is 0 if the method was not successful or 1 if it was
% successful, and iterations, which is the number of times the method was 
% iterated.

format long
syms x real
success = false;
count = 1;
r(1,1) = (b-a)*(f(a)+f(b))/2;
h = (b-a);
j = 1;

while (count < max && ~success)
   count = count + 1;
   j = j + 1;
   h = h / 2;
   fsum = 0;
   for i=1:2^(j-2)
       fsum = fsum + f(a + (2*i-1)*h);
   end;
   r(j,1) = r(j-1,1)/2 + h*fsum;
   for k=2:j
       r(j,k) = (4^(k-1)*r(j,k-1)-r(j-1,k-1))/(4^(k-1)-1);
   end;
   success = abs(r(j,j)-r(j-1,j-1)) < tol;
end;

ans = r(count,count);
flag = success
iterations = count