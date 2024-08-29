function error_code=pco_triggerCamera(camera)
if ~ismember(camera.CameraMode,[17 33])
    warning(['You tried to send a software trigger and the camera ' ...
        'isn''t in software mode you dumb dumb!']);
   return; 
end
[error_code] = pfTRIGGER_CAMERA(camera.BoardHandle);
if(error_code~=0) 
    error(['Could not trigger camera. Error is ',int2str(error_code)]);
    return;
end
end


