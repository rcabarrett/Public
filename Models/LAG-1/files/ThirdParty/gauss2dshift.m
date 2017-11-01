function g = gauss2dshift(kSize,mu,sigma)

[range_x,range_y]=meshgrid(-kSize:kSize);

g = exp(-0.5 * ((range_x-mu(1)).^2 + (range_y-mu(2)).^2) / sigma^2);



