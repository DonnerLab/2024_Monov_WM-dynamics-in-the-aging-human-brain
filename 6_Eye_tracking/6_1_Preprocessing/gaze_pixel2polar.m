function [pa,ecc]=gaze_pixel2polar(x,y,res,w,dist)
% compute pixels per degree of visual angle based on experimental setup
% res = [X X];  % monitor resolution [x,y]
% w = X;   % width of monior (in cm)
% dist = X;   % distance of participant's eyes from screen (in cm)

o = tan(0.5*pi/180)*dist;
ppd = 2*o*res(1)/w;   % pixels per degree of visual angle

% if gaze coordinates x and y are recorded in pixel positions, convert these to degrees of visual angle relative to fixation
x = (x-res(1)./2)./ppd;
y = (y-res(2)./2)./ppd;

% calculate polar angle and eccentricity of gaze position
[theta,r] = cart2pol(x,y);   % theta = polar angle, r = eccentricity
pa = theta; 
ecc = r; 
end 