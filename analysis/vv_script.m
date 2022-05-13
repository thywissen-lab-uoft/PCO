

%100Er spec before calibration
% paths = [ 
%         "Y:\Data\2022\2022.02\02.11\13 RF 75 spec 200Er 198G"
%         "Y:\Data\2022\2022.02\02.11\14 RF 75 spec 100Er 198G"
%         "Y:\Data\2022\2022.02\02.11\16 RF 75 spec 100Er 197G"
%         "Y:\Data\2022\2022.02\02.11\17  RF 75 spec 100Er 198G"
%         "Y:\Data\2022\2022.02\02.11\18  RF 75 spec 100Er 199G"
%         "Y:\Data\2022\2022.02\02.11\19  RF 75 spec 100Er 200G"
%         "Y:\Data\2022\2022.02\02.11\20  RF 75 spec 100Er 201G"
%         "Y:\Data\2022\2022.02\02.11\21  RF 75 spec 100Er  201.5G"
%         "Y:\Data\2022\2022.02\02.11\22  RF 75 spec 100Er 205G spin juggle loading"
%         "Y:\Data\2022\2022.02\02.11\23  RF 75 spec 100Er 206G spin juggle loading"
%         "Y:\Data\2022\2022.02\02.11\24  RF 75 spec 100Er 206G Ramp through resonance"
%         "Y:\Data\2022\2022.02\02.11\25  RF 75 spec 100Er 207G spin juggle loading"
%         "Y:\Data\2022\2022.02\02.11\26  RF 75 spec 100Er 208G spin juggle loading"
%         "Y:\Data\2022\2022.02\02.11\27 RAR"
%         ];
%       



%200Er spec attractive side after calibration
% paths = [ 
%         "Y:\Data\2022\2022.02\02.22\12 RF spec 200Er 203.5G"
%         "Y:\Data\2022\2022.02\02.22\13 RF spec 200Er 204G"
%         "Y:\Data\2022\2022.02\02.22\14 RF spec 200Er 204.5 G"
%         "Y:\Data\2022\2022.02\02.22\15 RF spec 200Er 208 G"
%         "Y:\Data\2022\2022.02\02.22\16 RF spec 200Er 209G"
%         ];


%200Er spec repuslive side and RAR
paths = [ 
%         "Y:\Data\2022\2022.02\02.18\16 95 spec 200Er 197G"
%         "Y:\Data\2022\2022.02\02.18\17 95 spec 200Er 198G"
%         "Y:\Data\2022\2022.02\02.18\18 95 spec 200Er 199G"
        "Y:\Data\2022\2022.02\02.18\19 95 spec 200Er 200G"
        "Y:\Data\2022\2022.02\02.18\20 95 spec 200Er 206G 4dBm SRS pow ramp through resonance"
        "Y:\Data\2022\2022.02\02.18\21 95 spec 200Er 206G 5dBm SRS pow ramp through resonance"
        "Y:\Data\2022\2022.02\02.18\22 95 spec 200Er 207G ramp through resonance"
        "Y:\Data\2022\2022.02\02.18\23 95 spec 200Er 208G ramp through resonance"
        ];
      

for vv=1:length(paths)
    imgdir = convertStringsToChars(paths(vv));
    pco_main
    
end
     
      


