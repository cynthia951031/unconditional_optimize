%problem 1
clc; clear; warning off;
f = @(x) (log(1 + x(1)^2) - x(2));
ce = @(x) (1 + x(1)^2)^2 + x(2) ^2 - 4;
ci = @(x) [];
xStar = [0, sqrt(3)];
[x, fval] = penalty_function_method(f, ce, ci, [2, 2]);
%[x, fval] = augmented_lagrangian_method(f, ce, ci, [2, 2]);
fprintf(1, 'problem 1: \n');
fprintf(1, 'final result: \n');
x
fval
fprintf(1, 'infinity norm number: \n');
normNum = norm(xStar - x, inf);
normNum

%problem 2
clc; clear;
f = @(x) ((x(1) - 1)^2 + (x(1) - x(2)) ^2 + (x(2) - x(3))^4);
ce = @(x) [x(1) * (1 + x(2)^2) + x(3)^4 - 4 - 3 * sqrt(2)]; 
ci = @(x) [x(1) + 10, 10 - x(1),...
    x(2) + 10, 10 - x(2), ...
    x(3) + 10, 10 - x(3)];
xStar = [1.104859024, 1.196674194, 1.535262257];
[x, fval] = penalty_function_method(f, ce, ci, [2, 2, 2]);
%[x, fval] = augmented_lagrangian_method(f, ce, ci, [2, 2, 2]);
fprintf(1, 'problem 2: \n');
fprintf(1, 'final result: \n');
x
fval
fprintf(1, 'infinity norm number: \n');
normNum = norm(xStar - x, inf);
normNum



%problem 3
clc; clear;
f = @(x) (x(1) * x(4) * (x(1) + x(2) + x(3)) + x(3));
ce = @(x) [x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2 - 40];
ci = @(x) [x(1)*x(2)*x(3)*x(4) - 25, ...
    x(1) - 1, 5 - x(1), ...
    x(2) - 1, 5 - x(2), ...
    x(3) - 1, 5 - x(3), ...
    x(4) - 1, 5 - x(4)];
xStar = [1, 4.7429994, 3.8211503, 1.3794082];
%[x, fval] = penalty_function_method(f, ce, ci, [1, 5, 5, 1]);
[x, fval] = augmented_lagrangian_method(f, ce, ci, [1, 5, 5, 1]);
fprintf(1, 'problem 3: \n');
fprintf(1, 'final result: \n');
x
fval
fprintf(1, 'infinity norm number: \n');
normNum = norm(xStar - x, inf);
normNum


