function passTolCheck = checkWithinTolerance(actual, expected, tolerance, checkInIfTrue)

passTolCheck = false;
% flag to see if the one must stay within the tolerance or be outside
% the tolerance
checkOutIfTrue = ~checkInIfTrue;
if checkInIfTrue && (abs(actual - expected) < tolerance)
    passTolCheck = true;
elseif checkOutIfTrue &&  (abs(actual-expected) > tolerance)
    passTolCheck = true;
end
