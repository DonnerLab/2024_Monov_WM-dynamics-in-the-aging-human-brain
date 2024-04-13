% Setup various options
function options = setup
% options.et = 'no';                % use EyeLink?
% options.do_trigger = 'no';        % send triggers via parallel port?
% options.kb_setup = 'psychophys';  % response setup: if 'MEG', will register multiple different input devices
% 
% options.dist = 60; % viewing distance in cm 
% options.width = 34.5; % physical width of the screen in cm
% options.height = 19; % physical height of the screen in cm
% 
% options.resolution = [1920, 1080];
% options.ppd = estimate_pixels_per_degree(options);

%%%%%%%%%%%% MEG SETTINGS
options.et = 'yes';                % use EyeLink?
options.do_trigger = 'yes';        % send triggers via parallel port?
options.kb_setup = 'MEG';  % response setup: if 'MEG', will register multiple different input devices

options.dist = 61; % viewing distance in cm 
options.width = 51.9347; % physical width of the screen in cm
options.height = 29.5; % physical height of the screen in cm

options.resolution = [1920, 1080];
options.ppd = estimate_pixels_per_degree(options);
end

function ppd = estimate_pixels_per_degree(options)
o = tan(0.5*pi/180)*options.dist;
ppd = 2*o*options.resolution(1)/options.width;
end