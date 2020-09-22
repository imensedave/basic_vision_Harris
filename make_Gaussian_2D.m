function mask = make_Gaussian_2D(sigma)
% function to make a 2D Gaussian mask 
% copyright d.sinclair 2020

nx = ceil( 2.2 * sigma +0.5 );
num = 2*nx + 1;

mask = zeros(num,num);
sigma2 = 2*sigma*sigma;

nx = nx+1;

for i = 1: nx
    X = i-nx;
    for j = 1: nx
        Y = j-nx; 
        globit = 0;
        nu = 0;
        for x=-0.5:0.1:0.5
            xx = X+x;
            for y=-0.5:0.1:0.5
                yy = Y+y;
                globit = globit + exp( - (xx*xx + yy*yy)/ (sigma2) ) ;
            end
        end
        mask(i,j) = globit;
    end
end

m2 = flipdim(mask,1);
mask = max(mask,m2);
m3 = flipdim(mask,2);
mask = max(mask,m3);

s = sum(mask(:));
mask = mask/s;
   
