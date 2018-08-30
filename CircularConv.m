function x_ = CircularConv(x, ksize, poly_is_closed)
%{
-About-
convolves a vector circularly with a specified kernel size. Circularly
means that the convolution wraps around the vector from the end to the
begining. 

-Inputs-
x:  1-dimensional vector to be convolved
ksize: size of the convolving kernel, must be an odd integer. If ksize is
not an odd integer it will automatically (and silently) be converted to the
next sequential odd integer.
poly_is_closed: a bool specifying whether or not the polygon is closed

-varargin-
'error_lim'     set the limit of acceptable error in the analysis

-Outputs-
output:     result of x,y projected into some new space

-Example-
x = [1:10, 10:-1:1, 1];
y = zeros(size(x));
y(1:10) = 5;
y(end) = 5;

plot(x,y)
plot(CircularConv(x,7,true), CircularConv(y,7,true), 'r')
   
-Supplementary-

-Keywords-
vector polygon convolution

-Dependencies-

-References-

-Author-
Brad Parry, 2018 May 24
%}



ksize = floor(ksize/2)*2+1;

x=x(:)';
dy = floor(ksize/2);

if poly_is_closed
    x = x(1:end-1);
end

x_ = [x(length(x)-(dy-1):end),x, x(1:dy)];

x_ = conv(x_,ones(ksize,1)/ksize,'valid');

if poly_is_closed
    x_(end+1) = x_(1);
end
