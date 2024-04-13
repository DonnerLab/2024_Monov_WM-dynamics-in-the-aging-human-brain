% displays brief instructions screen at the start of each block, specifying
% current block number

function block_instructions(b,window,windowRect,language)

% Get window center
[~,Y] = RectCenter(windowRect);

% Specify text
if strcmp(language,'german')
    str1 = sprintf('Sie starten jetzt mit Block Nummer %s.',b);
    str2 = 'Drücken Sie einen Knopf zum Starten...';
elseif strcmp(language,'english')
    str1 = sprintf('You will now start block number %s.',b);
    str2 = 'Press a button to begin...';
end

% Display text
DrawFormattedText(window,str1,'center',Y-40,[255 255 255]);
DrawFormattedText(window,str2,'center',Y+40,[255 255 255]);

Screen('Flip', window);

% Wait for response
WaitSecs(0.5);
keyIsDown = 0;
while ~keyIsDown
    [keyIsDown, ~, ~] = KbCheck;  % checking for response
end