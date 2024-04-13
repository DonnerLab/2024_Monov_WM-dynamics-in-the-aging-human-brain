function end_block_screen(block,loadstr,window,windowRect,language)

% Calculate current block and total earnings
acc = 0; ntrials = 0;
for b = 1:block;
    load([loadstr,num2str(b),'.mat']);
    acc = acc+nansum(Behav(:,6));
    ntrials = ntrials+size(Behav,1);
end

% Get window center
[~,Y] = RectCenter(windowRect);

% Specify text
if strcmp(language,'german')
    str1 = 'Block beendet!';
    str2 = sprintf('Richtige Antworten in diesem Block = %2.1f Prozent',nansum(Behav(:,6))/size(Behav,1)*100);
    str3 = sprintf('Richtige Antworten heute = %2.1f Prozent',acc/ntrials*100);
    str4 = sprintf('Sie haben %d mal die Fixierung verloren.',length(find(Behav(:,end)>0)));
    str5 = 'Versuchen Sie bitte ausschließlich während der Blinzel Periode ';
    str6 = 'zwischen den Aufgaben das Fixationskreuz nicht zu fixieren!';
    str7 = 'Bitte warten Sie während wir Ihre Daten speichern...';
elseif strcmp(language,'english')
    str1 = 'Block complete!';
    str2 = sprintf('Accuracy for that block = %2.1f percent',nansum(Behav(:,6))/size(Behav,1)*100);
    str3 = sprintf('Total accuracy for this session = %2.1f percent',acc/ntrials*100);
    str4 = sprintf('You broke fixation during a trial on %d trials.',length(find(Behav(:,end)>0)));
    str5 = 'Remember to try as hard as you can to only do this during the blink period between trials!';
    str6 = '';
    str7 = 'Please wait while we save your data...';
end

% Display text
DrawFormattedText(window,str1,'center',Y-250,[255 255 255]);
DrawFormattedText(window,str2,'center',Y-110,[255 255 255]);
DrawFormattedText(window,str3,'center',Y-60,[255 255 255]);
DrawFormattedText(window,str4,'center',Y+80,[255 255 255]);
DrawFormattedText(window,str5,'center',Y+130,[255 255 255]);
DrawFormattedText(window,str6,'center',Y+180,[255 255 255]);
DrawFormattedText(window,str7,'center',Y+300,[255 255 255]);

Screen('Flip', window);

WaitSecs(5);