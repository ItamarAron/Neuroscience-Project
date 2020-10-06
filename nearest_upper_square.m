function [y] = nearest_upper_square(x)
%Nearest upper square takes a number x and finds the closest perfect square that is greater than x.
while mod(x,sqrt(x)) ~= 0
    x = x+1;
end
y=x;
end

