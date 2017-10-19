function G = gauss2d_LAG-1(half_range, sigma, varargin)
  
    half_range_x = half_range;
    half_range_y = half_range;
    mu_x = 0;
    mu_y = 0;

    if ~isempty(varargin)
        half_range_x = varargin{1};
        half_range_y = varargin{2};
        mu_x = varargin{3};
        mu_y = varargin{4};
    end

        [x,y] = meshgrid(-half_range_x:half_range_x,-half_range_y:half_range_y);

        G = exp( -0.5*((x-mu_x).^2 + (y-mu_y).^2)/sigma^2);
end


%{
function g = gauss2d(kSize, sigma)

[range_x,range_y]=meshgrid(-kSize:kSize);

g = exp(-0.5 * (range_x.^2 + range_y.^2) / sigma^2);
%}