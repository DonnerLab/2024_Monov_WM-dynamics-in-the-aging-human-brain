function [imgout1] = make_gabor_mask_tex(ppd,res,g_r,g_e)

% ppd      =  pixels per d.v.a
% res      =  resolution of display [horizontal, vertical]
% g_r      =  radius (d.v.a.) of gabors
% g_e      =  eccentricity from fixation of memoranda

dva_res = res./2./ppd;   % screen resolution expressed as d.v.a.
g_xy = sqrt((g_e.^2)/2); % converting eccentricity into horizontal/vertical offset

% make transparent mask for 4 possible eccentric gabors (i.e. memoranda)
[x1, y1] = meshgrid(linspace(-dva_res(1)+g_xy, dva_res(1)+g_xy, res(1)), linspace(-dva_res(2)+g_xy, dva_res(2)+g_xy, res(2)));  % create matrix of X and Y d.v.a. values for each pixel RELATIVE TO TOP-LEFT PATCH
r1 = sqrt(y1.^2 + x1.^2);  % calculate polar radius coordinate for each pixel RELATIVE TO LEFT PATCH
mask1 = r1>g_r;

[x2, y2] = meshgrid(linspace(-dva_res(1)-g_xy, dva_res(1)-g_xy, res(1)), linspace(-dva_res(2)+g_xy, dva_res(2)+g_xy, res(2)));  % create matrix of X and Y d.v.a. values for each pixel RELATIVE TO TOP-RIGHT PATCH
r2 = sqrt(y2.^2 + x2.^2);  % calculate polar radius coordinate for each pixel RELATIVE TO RIGHT PATCH
mask2 = r2>g_r;

[x3, y3] = meshgrid(linspace(-dva_res(1)-g_xy, dva_res(1)-g_xy, res(1)), linspace(-dva_res(2)-g_xy, dva_res(2)-g_xy, res(2)));  % create matrix of X and Y d.v.a. values for each pixel RELATIVE TO BOTTOM-RIGHT PATCH
r3 = sqrt(y3.^2 + x3.^2);  % calculate polar radius coordinate for each pixel RELATIVE TO RIGHT PATCH
mask3 = r3>g_r;

[x4, y4] = meshgrid(linspace(-dva_res(1)+g_xy, dva_res(1)+g_xy, res(1)), linspace(-dva_res(2)-g_xy, dva_res(2)-g_xy, res(2)));  % create matrix of X and Y d.v.a. values for each pixel RELATIVE TO BOTTOM-LEFT PATCH
r4 = sqrt(y4.^2 + x4.^2);  % calculate polar radius coordinate for each pixel RELATIVE TO LEFT PATCH
mask4 = r4>g_r;

unigrey = ones(size(r1)).*127;  % initializing pure grey image
alphaout = ones(size(r1)).*255;   % intializing matrix of transparency values
alphaout(mask1==0 | mask2==0 | mask3==0 | mask4==0) = 0;  % setting everything within radius of each possible gabor to be completely transparent (0 = fully transparent, 255 = full opaque)
imgout1 = cat(3,unigrey,unigrey,unigrey,alphaout);

