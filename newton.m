function [ root, flag, iterations ] = newton(f, df, x0, max, tol)
%f is the input function, df is the derivative of the input function
%x0 is the initial guess, and max is the maximum number of iterations
%to perform. tol is the relative tolerance to test whether the method
%has reached success.
%
%The function returns three outputs: root, the root if the test was
%successful, or the last guess if the test was not successful, flag,
%which is 0 if the test was a failure and 1 if it was a success, and
%iterations, which is simply the number of iterations executed during
%the test.
format long
syms x real
success = false;
count = 0;
 while (count < max) && ~success
     x1 = x0;
     x0 = x0 - f(x0)/df(x0);
     count = count + 1;
    if abs((x0-x1)/x0) < tol
          success = true;
     end;
 end;
 if success
     root = x1
     flag = true
     iterations = count - 1
 else
     root = x0
     flag = false
     iterations = count
 end