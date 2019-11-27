function TF = gazeIsInWindow(relativeEyePosition,r)
%gazeIsInWindow: a function for Ex to determine if the gaze is in a
%specific window. -Adam C. Snyder adam@adamcsnyder.com
%
%   relativeEyePosition should be the eye position relative to the center
%   of the desired window (e.g., eyePos-[fixX fixY] in the calling
%   function).
%
%   r should be a scalar or a 2x1 vector. If it is a scalar, a circular
%   window is checked. If it is a 2x1 vector, a rectangular window is
%   checked. Square windows are just a special case of a rectangular
%   window.
%
%   The output is logical true if the gaze is in the window and logical
%   false if the gaze is outside the window or exactly on the window
%   border.
%
%   NB: I've changed waitForFixation so that it no longer calls this
%   function -acs09dec2015

    switch size(r,1)
        case 1 %indicates a circular window is desired
            TF = sum(relativeEyePosition.^2)<r.^2;
        case 2 %indicates a rectangular window is desired
            TF = all(abs(relativeEyePosition(:))<abs(r(:)));
        otherwise
            error('Unsupported window shape: ''r'' is expected to be a 1xNwindows (circular) or 2xNwindows (rectangular) array');
    end;        
        
