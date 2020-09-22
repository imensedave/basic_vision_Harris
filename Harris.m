function [pts] = Harris(I,sigma, k)
%
%
%
%
%   Harris corner fudge....
%   cim = imread('image.jpg');
%   I = double(rgb2gray(cim(1:400,1:400,:)));
%   pts = Harris(I, 2.0, 0.04);
%
% copyright d.sinclair 2020


%sigma=2; 
dx = [-1 0 1; -1 0 1; -1 0 1]; % The gradient operator
dy = dx';

Ix = conv2(I, dx, 'same');
Iy = conv2(I, dy, 'same');
g = make_Gaussian_2D(sigma);


Ix2 = conv2(Ix.^2, g, 'same');
Iy2 = conv2(Iy.^2, g, 'same');
Ixy = conv2(Ix.*Iy, g,'same');

%k = 0.04;
%Hassisiness:
H = (Ix2.*Iy2 - Ixy.^2) - k*(Ix2 + Iy2).^2;

d = ceil(size(g,1)/2) + 3;


% non max suppression

[nr,nc] = size(I);
br = nr-d;
bc = nc-d;
ncx = nc-2*d+1;

pts = zeros(0,3);
cnt = 0;

for r=d:nr-d
    R0 = H(r-1,d:bc);
    R = H(r,d:bc);
    R1 = H(r+1,d:bc);
    ids = find( R(2:ncx-1) > 200 & R(2:ncx-1) > R0(1:ncx-2) & R(2:ncx-1) > R0(3:ncx) & R(2:ncx-1) > R0(2:ncx-1) &... 
    R(2:ncx-1) > R(1:ncx-2) & R(2:ncx-1) >= R(3:ncx)  &...
    R(2:ncx-1) >= R1(1:ncx-2) & R(2:ncx-1) >= R1(3:ncx) & R(2:ncx-1) >= R1(2:ncx-1)       );
    
    cs = ids-1+d;
    nmx = size(cs,2);
    for z=1:nmx
        cnt = cnt+1;
        pts(cnt,1) = r;
        pts(cnt,2) = cs(z);
        pts(cnt,3) = H(r,cs(z));
    end
end

[vals, ord] = sort(pts(:,3),'descend');

pts = pts(ord(:),:);

nump = size(pts,1)

figure( 8)
imagesc( I)
hold on
if( nump > 2000 )
    nump = 2000;
end
plot( pts(1:nump,2), pts(1:nump,1), 'rx')
hold off

