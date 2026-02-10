This repository contains four MatLab programs to run approxiamtion methods for function roots in differential equations.
The four methods: Bisection, Fixed-Point-Iteration, Secant, and Newton's Method were all carried out on the example function:
f(x)=(1/3)*sin(x) - x - (1/5), with the exmaple codes for each shown below.

## Example: Bisection Method

```Matlab

% Bisection Method
% Stopping criterion: |x_k - x_{k-1}| < 1e-8

f = @(x) (1/3)*sin(x) - x - 1/5;   % define function f
a = -0.5;     % assign A with left endpoint
b = 0;        % assign b with right endpoint
tol = 1e-8;   % assign t with tolerance
k = 0;   %initialize variable k which tracks the number of iterations
x_old = (a + b)/2; %create first bisection estimate
while true %create while loop to make estimates
    k = k + 1; %iterate k 
    if f(a)*f(x_old) < 0
        b = x_old; %assign new b value
    else
        a = x_old;%assign new a value
    end
    x_new = (a + b)/2; %get new estimate
    fprintf('%3d   %.9f\n', k, x_new);% Print iteration number and estimate
    if abs(x_new - x_old) < tol %assign stopping criteria
        break;
    end
    x_old = x_new; %replace the old value so we can restart the while loop
end

```

## Fixed-Point-Iteration Method
```Matlab

% Fixed Point Iteration Method
% Stopping criterion: |x_k - x_{k-1}| < 1e-8

g = @(x) (1/3)*sin(x) - 1/5; %create g as our example fixed-point function
x0 = 0;    % set x0 to our initial estimate
newtol = 1e-8;   % set newtol = the stopping condition
j = 0; %create j to count the iterations
xold = x0; %set the first old x to x0
while true %create iteration while loop
    j = j + 1; %iterate j
    xnew = g(xold);%set new x value to the result of g(xold)
    fprintf('%3d   %.9f\n', j, xnew);%print iteration number and estimate
    if abs(xnew - xold) < newtol %create stopping condition
        break;
    end
    xold = xnew;
end
```

## Newton's Method
```Matlab
%Newton's method

function [x] = newton(f, df, x0, tol) %create newton method function
    x = x0; % Initialize x with x0 initial guess
    xold = x0; % Initialize xold with x0
    k = 0;
    while true
        k = k + 1;
        fx = f(x);
        f1x = df(x);
        xnew = x - (fx/f1x);
        if abs(xnew - x) < tol
            x = xnew;
            return;
        end
        x = xnew;
        % Print iteration number and iterate
        fprintf('%3d   %.9f\n', k, xnew);
        % Stopping criterion
        if abs(xnew - xold) < tol
            break;
        end
        xold = xnew;
    end
end

%example usage on the same function f(x) = (1/3)sin(x)-x-(1/5)

f = @(x) (1/3)*sin(x) - x - 1/5;  % define function f
df = @(x) (1/3)*cos(x)-1; % define df as derivative of f
df2 = @(x) (-1*1/3)*sin(x);
tol = 1e-10;   % assign tol with tolerance
k = 0;   %initialize variable k which tracks the number of iterations
x0 = 0; %create first guess for x0

newton(f, df, x0, tol) %run function 

```

## Secant Method
```Matlab

% Secant Approximation Method
% Stopping criterion: |x_k - x_{k-1}| < 1e-10
function x = secant_method(f, x0, x1, tol) % create secant method function
   k = 0;
   while true
       k = k + 1;  % Iterate k
       fx0 = f(x0); %assign fx0 with f(x0) and fx1 with f(x1)
       fx1 = f(x1);
       xnew = x1 - ((fx1*(x1 - x0))/(fx1 - fx0)); % xnew becomes next xk
       fprintf('%3d   %.10f\n', k, xnew); % print k and xk with up to 10 decimals
       if abs(xnew - x1) < tol %create stopping criteria
           x = xnew;
           break;
       end
       x0 = x1; %assign new x0 and x1 values to restart the loop
       x1 = xnew;
   end
end
f = @(x) (1/3)*sin(x) - x - 1/5;  % define function f
tol = 1e-10;   % assign tol with tolerance level
x0 = 0;        % assign x0 with lower first estimate 0
x1 = 1;    %assign x1 with higher first estimate 1
secant_method(f, x0, x1, tol)
```
## Conclusion
Newton's method converged the fastest to the true root x = -0.2978. This is because Newton's method converges quadrilatically, whereas the secant method converges superlinearly, and bisection and fixed-point-iteration converge linearly.
