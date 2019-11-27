function shouldBreak = keyboardEvents()
global trialMessage;

shouldBreak = 0;

[ keyIsDown, keyCode] = KbQueueCheck;
if keyIsDown
    c = KbName(keyCode);
    KbQueueFlush;
    if numel(c)>1, return; end;
    %    if CharAvail
    %        c = GetChar;
    switch c
        case 'j'
            sendCode(18);
            giveJuice();
        case 'q'
            trialMessage = -1;
            shouldBreak = 1;
        case 'E'
            error('keyboardEvents:errorRequested','Error requested at keyboard. This is intended to be used for testing purposes.');
    end
end

end
