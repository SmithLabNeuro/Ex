% diode_loop.m

while CharAvail;
    GetChar;
end

%% screen initialization
Screen('Preference', 'SkipSyncTests', 0);

%gam = load(photometerFile);
%vals = makeGammaTable(gam(:,1),gam(:,2));
%creen('LoadNormalizedGammaTable', 0, vals);

%% UDP initialization
delete(instrfind)
%u = udp(ipAddress,'RemotePort',8844,'LocalPort',8866);
%fopen(u);
%local_ip = char(java.net.InetAddress.getLocalHost.getHostAddress);
matlabUDP('close');
matlabUDP('open','192.168.1.10', '192.168.1.11', 4243);

%% Other initialization
% specifies size and location of photodiode square:
% upper left of screen, 50 x 50 pixels
diodeLoc = [10 10 40 40]; %made this slightly smaller to prevent light leaking out the side. -ACS 17Jul2013

% bottom right of screen, 65 x 65 pixels
%diodeLoc = [960 704 1024 768];

AssertOpenGL;

screens = Screen('Screens');
screenNumber = max(screens);

sv.white = WhiteIndex(screenNumber);
sv.black = BlackIndex(screenNumber);
sv.gray = (sv.white+sv.black)/2;
if round(sv.gray) == sv.white;
    sv.gray = sv.black;
end
sv.inc = sv.white-sv.gray;
sv.bgColor = sv.gray;
diodeVal = 0;

[w sv.screenRect] = Screen('OpenWindow',screenNumber,sv.bgColor);
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

sv.midScreen = sv.screenRect(3:4)/2;

% show the photodiode square
Screen(w,'FillRect',sv.black,diodeLoc);
Screen('Flip',w);

% turn off mouse pointer
HideCursor();

% get the frame time (inter-frame interval)
%sv.ifi = Screen('GetFlipInterval',w);

%% DISPLAY LOOP
while(1) % display loop

    if matlabUDP('check') % received network message
        [s1 s] = strtok(matlabUDP('receive'));
        
        %for debugging use this line
        disp([s1 s]);
        
        switch s1
            case 'diodeon'
                % arguments: object ID to tie the diode square to
                diodeVal = 255;
            case 'diodeoff'
                % arguments: object ID to tie the diode square to
                diodeVal = 0;
        end
        
        matlabUDP('send', s1);
    end
        
    % display section
    Screen(w,'FillRect',diodeVal * [1 1 1],diodeLoc);
    Screen('Flip',w);
    
    % check for keyboard input
    if CharAvail
        c = GetChar;
        if c == 'x'
            break;
        end
    end
end;

sca
