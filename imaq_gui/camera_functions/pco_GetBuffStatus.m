
function buff_status = pco_GetBuffStatus(camera,n)
% Checks the status of a buffer
%
%   camera     - the primary handles object
%   n          - which buffer to check (1,2, or 3)
%   buff_status - the status of the buffer
%                   1 queued
%                   2 not queued
%                   3 waiting with image
%                   4 getting image 
%                   -1 not allocated

% Assume the buffer is not allocated
buff_status = -1;

if camera.buf_nums(n)~=-1    
    % Read the buffer
    [error_code,buff_status] = pfGETBUFFER_STATUS(...
        camera.BoardHandle,camera.buf_nums(n),0,4);
        
    % Process the error code if read fails
    if error_code 
        error_code = pco_uint32err(error_code);            
        error(['Could not determine buffer status. Error is ' ...
            '0x' dec2hex(error_code)]);
     	return;
    end
    
    % Convert it into readable value
    buff_status=pco_uint32(buff_status);    

    % Look at bits associated with the parts we want
    if (bitand(buff_status,4))
        buff_status = 1;
    elseif (bitand(buff_status,2))
        buff_status = 3;
    end    
    % Display the buffer status (for debugging)
%       disp(['Buffer (' num2str(n) ') : ' ...
%           num2str(camera.buf_nums(n)) ' status is ',num2str(buff_status,'%08X')]);
end
    
end


