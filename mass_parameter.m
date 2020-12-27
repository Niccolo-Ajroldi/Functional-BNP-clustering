clc
clear all


syms a
expr = 1;
for i = 2:25
    expr = expr + a/(a+i);
end

vpasolve(expr-2,a)

vpasolve(a*log(1+26/a)-2,a)
