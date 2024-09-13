
function camera=pco_initCam(camera)   
    fprintf('Intializing PCI PCO 540 board and camera ...');

    % If no arguments, go to default settings    
    if nargin==0
        camera=initCamStruct;
    end
    
    % Reset camera
    camera.RunStatus=0;         % Camera status (1 = on, 0 = off);
    camera.NumAcquired=0;       % Number of acquired images (0,1,2)
    
    % Connect to the PCI PCO 540 Board
    [error_code,camera.BoardHandle] = pfINITBOARD(0);    
    if  error_code~=0
        warning('Unable to initialize board!!');
        error(['Could not connect to the board. Error is ' ...
            '0x' dec2hex(pco_uint32err(error_code))]);
        return;
    end       
    camera.isConnected=1;
    [error_code, value] = pfGETBOARDVAL(camera.BoardHandle,'PCC_VAL_BOARD_STATUS');
    if(error_code)
        pco_errdisp('pfGETBOARDVAL',error_code);    
    else
        if(bitand(value,hex2dec('01'))==hex2dec('01'))            
            disp('Camera is running call STOP_CAMERA')     
            error_code=pfSTOP_CAMERA(camera.BoardHandle);
            pco_errdisp('pfSTOP_CAMERA',error_code);
        end 
    end        
    disp('Done');    
    % Configure the camera
    camera=configCam(camera);    
end