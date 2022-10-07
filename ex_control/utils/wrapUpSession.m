function wrapUpSession()

global sessionNumber sqlDb params

if ~isempty(sqlDb)
    subject = params.SubjectID;
    waterAndThresh = inputdlg({'Milliliters water for session:','Trellis RMS spike threshold for session:'});
    
    updateExperimentSessionInDatabase(sqlDb, sessionNumber, subject, 'water_ml', waterAndThresh{1}, 'threshold_rms', waterAndThresh{2});
end

socketsDatComp.sender = matlabUDP2('open',params.control2dataIP,params.data2controlIP,params.control2dataSocketSend);
socketsDatComp.receiver = matlabUDP2('open',params.control2dataIP,params.data2controlIP,params.control2dataSocketReceive);

try
    rc = sendMessageWaitAck(socketsDatComp, 'sessionEnd', 1);
catch err
    if any(strfind(err.message, 'Communication with data computer failed--is it running recordFromControl.m?'))
        fprintf('%s\n', err.message);
    end
end
matlabUDP2('all_close');
end