function success = joystickGetPosition(~, ~, analogChannelPosInput)

posOut = zeros(size(analogChannelPosInput));
for channel = 1:length(analogChannelPosInput)
    chPos = unixGetAnalogInput(analogChannelPosInput);
    posOut(channel) = chPos;
end

success = posOut;