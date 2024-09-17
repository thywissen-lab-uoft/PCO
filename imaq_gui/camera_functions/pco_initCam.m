
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

%% Define configCam function


function camera=configCam(camera)
    disp(' ');
    disp('Configuring camera acquisition');
    disp(['     ExposureTime : ' num2str(camera.ExposureTime) ' us']);
    disp(['     CameraMode   : ' num2str(camera.CameraMode)]);
    disp(['     NumImages    : ' num2str(camera.NumImages)]);

    % Set the camera mode (mode, exposure, binning, etc.)
    bit_pix=12;
    auto_exp=50;    % auto exposure level (ignored)
    [error_code] = pfSETMODE(camera.BoardHandle,...
        camera.CameraMode,...
        auto_exp,...
        camera.ExposureTime,...
        0,0,0,0,bit_pix,0);      

    if (error_code)
        pco_errdisp('pfSETMODE',error_code)
        warning('Oh no something went wrong on configuring BAD.'); 
    end
    

    % Read in the camera image size
    [error_code,ccd_width,ccd_height,act_width,act_height,bit_pix]=...
        pfGETSIZES(camera.BoardHandle);
    if (error_code)
        pco_errdisp('pfGETSIZES',error_code)
        warning('oh no something went wrong on reading the sensor.'); 
    end

    camera.H=double(act_height);
    camera.W=double(act_width);
    camera.BitDepth=bit_pix;    

    % Determine size of buffer to allocate 
    imasize=act_width*act_height*floor((bit_pix+7)/8);  
    image_stack=ones(act_width,act_height,camera.NumImages,'uint16');      

    % Allocate and map buffers for each image
    fprintf(['Allocating (' num2str(camera.NumImages) ') buffers ... ']);    
    for i = 1:camera.NumImages   
        fprintf([num2str(i) ' ']);
        ima_ptr(i) = libpointer('uint16Ptr',image_stack(:,:,i));
        [error_code, buf_nums(i)] = pfALLOCATE_BUFFER_EX(...
            camera.BoardHandle,-1,imasize,ima_ptr(i)); 
        fprintf(['(bnum=' num2str(buf_nums(i)) ') ']);        
    end
    disp('done');
    
    camera.buf_ptrs=ima_ptr;
    camera.buf_nums=buf_nums;
    
    % Remove buffers from list write (ie. clear them)
    fprintf('Clearing buffer queue ...');
    [error_code] = pfREMOVE_ALL_BUFFER_FROM_LIST(camera.BoardHandle);
    camera.NumAcquired=0;
    disp(' done');   
    
    % Read the anticpated CCD readout time (useful for double shutter
    % exposure)
    [error_code, rTime] = pfGETBOARDVAL(camera.BoardHandle, 'PCC_VAL_READOUTTIME');
    if ~error_code
        disp(['Calculated CCD read time per image is ' ...
            num2str(1E-3*double(rTime)) ' ms.']);
    end
    
    disp('Camera acquistion configured.');
    
    
end