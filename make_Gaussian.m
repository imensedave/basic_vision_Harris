function mask = make_Gaussian(sigma)
% function to make a 1D Gaussian mask 
%keyboard

taps = floor( 2.2 * sigma +0.5 );
mask = zeros(1,2*taps +1 );
sigma2 = 2*sigma*sigma;

for i = 0: taps
   tpx = i-0.44;
   globit = 0;
	for j=0:8
	  globit = globit + exp( - tpx * tpx/ (sigma2) ) ;
	  tpx = tpx + 0.11;
   end
   mask(taps+1-i) = globit;
   mask(i+taps+1) = globit;
   
end

s = sum(mask);
mask = mask/s;   
   