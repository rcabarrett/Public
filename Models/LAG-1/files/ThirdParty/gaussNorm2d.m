function G = gaussNorm2d(half_range, sigma, varargin)
  
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

    if any(G) %if any value is above 0
      G = G / sum(G(:));
    end
end