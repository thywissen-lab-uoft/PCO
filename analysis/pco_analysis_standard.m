%% pco_analysis_standard
%
% This script runs the basic plotting and analysis code for box count,
% gaussian, and erf fits.

%% Box Count Analysis
% This is the default box count analysis.

if doBoxCount
    boxPopts = struct;
    boxPopts.FigLabel = FigLabel;
    boxPopts.xUnit=pco_unit;
    boxPopts.NumberExpFit = 0;        % Fit exponential decay to atom number
    boxPopts.NumberLorentzianFit=0;   % Fit atom number to lorentzian
    boxPopts.CenterSineFit = 0;       % Fit sine fit to cloud center
    boxPopts.CenterDecaySineFit = 0;  % Fit decaying sine to cloud center
    boxPopts.CenterParabolaFit = 0;
    boxPopts.CenterLinearFit = 0;     % Linear fit to cloud center
    boxPopts.NumberExpOffsetFit = 0; % Exp decay fit with nonzero offset    
       
    hF_number_box = showAtomNumber(box_data,pco_xVar,boxPopts);  
    ylim([0 max(get(gca,'YLim'))]);    
    if doSave;saveFigure(hF_number_box,'box_number',saveOpts);end
    
    % Plot the ratios if there are more than one ROI.
    if size(ROI,1)>1    
        hF_number_box_ratio=showNumberRatio(box_data,pco_xVar,boxPopts);
        if doSave;saveFigure(hF_number_box_ratio,'box_number_ratio',saveOpts);end
    end

    % Cloud centre
    hF_Centre_box=showAtomCentre(box_data,pco_xVar,boxPopts);    
    if doSave;saveFigure(hF_Centre_box,'box_position',saveOpts);end 
        
    % box Size
    hF_size_box=showSize(box_data,pco_xVar,boxPopts);    
    if doSave;saveFigure(hF_size_box,'box_size',saveOpts);end       
end

%% 2D Gauss Analysis

if doGaussFit
% This is the default gaussian analysis.
    
    gaussPopts = struct;
    gaussPopts.FigLabel = FigLabel;
    gaussPopts.xUnit=pco_unit;
    gaussPopts.NumberExpFit = 0;        % Fit exponential decay to atom number
    gaussPopts.NumberLorentzianFit=0;   % Fit atom number to lorentzian
    gaussPopts.CenterSineFit = 0;       % Fit sine fit to cloud center
    gaussPopts.CenterDecaySineFit = 0;  % Fit decaying sine to cloud center
    gaussPopts.CenterParabolaFit = 0;
    gaussPopts.CenterLinearFit = 0;     % Linear fit to cloud center
    gaussPopts.NumberExpOffsetFit = 0; % Exp decay fit with nonzero offset
    
    % Plot the statistics of gaussian fit
    hF_stats=showGaussStats(gauss_data,gaussPopts);     
    if doSave;saveFigure(hF_stats,'gauss_stats',saveOpts);end       
    
    hF_number_gauss = showAtomNumber(gauss_data,pco_xVar,gaussPopts);  
    ylim([0 max(get(gca,'YLim'))]);    
    if doSave;saveFigure(hF_number_gauss,'gauss_number',saveOpts);end
    
    % Plot the ratios if there are more than one ROI.
    if size(ROI,1)>1    
        hF_number_gauss_ratio=showNumberRatio(gauss_data,pco_xVar,gaussPopts);
        if doSave;saveFigure(hF_number_gauss_ratio,'gauss_number_ratio',saveOpts);end
    end
    
    % Cloud Error
    hF_Error=showError(gauss_data,pco_xVar,gaussPopts);    
    if doSave;saveFigure(hF_Error,'gauss_error',saveOpts);end  
   
    % Cloud centre
    hF_Centre=showAtomCentre(gauss_data,pco_xVar,gaussPopts);    
    if doSave;saveFigure(hF_Centre,'gauss_position',saveOpts);end 
        
    % Gauss Size
    hF_size=showSize(gauss_data,pco_xVar,gaussPopts);    
    if doSave;saveFigure(hF_size,'gauss_size',saveOpts);end
    
    % Single shot temperature analysis
    [hF_tempsingle,Tdata]=showGaussSingleTemperature(gauss_data,pco_xVar,gaussPopts);    
    if doSave;saveFigure(hF_tempsingle,'gauss_tempsingle',saveOpts);end     
        
    % Aspect ratio
    hF_ratio=showAspectRatio(gauss_data,pco_xVar,gaussPopts);    
    if doSave;saveFigure(hF_ratio,'gauss_ratio',saveOpts);end   

    % Peak gaussian density
    hF_density=showDensity(gauss_data,pco_xVar,gaussPopts);    
    if doSave;saveFigure(hF_density,'gauss_density',saveOpts);end    

    % Gaussian Temperature Analysis
    if isequal(pco_xVar,'tof') && size(gauss_data.Natoms,2)>2
        [hF_temp,fitX,fitY]=computeGaussianTemperature(gauss_data,gaussPopts);
        if doSave;saveFigure(hF_temp,'gauss_temp',saveOpts);end    
    end      
    
    % Determine which ROIs to perform BEC analysis on (for double shutter)
    if doBEC 
        BECopts=struct;
        BECopts.FigLabel = FigLabel;
        BECopts.xUnit=pco_unit;   
        BECopts.pow2freq = @(P) 0.725*61.5*sqrt(P./(0.085)); % Calibrated 2021.02.25
        
        [hF_BEC,BECdata]=BECanalysis(gauss_data,pco_xVar,BECopts);    

        if doSave;saveFigure(hF_BEC,'gauss_BEC',saveOpts);end        
    end         
end

%% 2D Erf Analysis

if doErfFit  
% This is the default erf analysis.
    
    ErfPopts = struct;
    ErfPopts.FigLabel = FigLabel;
    ErfPopts.xUnit=pco_unit;
    ErfPopts.NumberExpFit = 0;        % Fit exponential decay to atom number
    ErfPopts.NumberLorentzianFit=0;   % Fit atom number to lorentzian
    ErfPopts.CenterSineFit = 0;       % Fit sine fit to cloud center
    ErfPopts.CenterDecaySineFit = 0;  % Fit decaying sine to cloud center
    ErfPopts.CenterParabolaFit = 0;
    ErfPopts.CenterLinearFit = 0;     % Linear fit to cloud center
    ErfPopts.NumberExpOffsetFit = 0; % Exp decay fit with nonzero offset
    
    % Plot the statistics of erfian fit
    hF_stats_erf=showErfStats(erf_data,ErfPopts);     
    if doSave;saveFigure(hF_stats_erf,'erf_stats',saveOpts);end       
    
    hF_number_erf = showAtomNumber(erf_data,pco_xVar,ErfPopts);  
    ylim([0 max(get(gca,'YLim'))]);    
    if doSave;saveFigure(hF_number_erf,'erf_number',saveOpts);end
    
    % Plot the ratios if there are more than one ROI.
    if size(ROI,1)>1    
        hF_number_erf_ratio=showNumberRatio(erf_data,pco_xVar,ErfPopts);
        if doSave;saveFigure(hF_number_erf_ratio,'erf_number_ratio',saveOpts);end
    end
    
    % Cloud Error
    hF_Error=showError(erf_data,pco_xVar,ErfPopts);    
    if doSave;saveFigure(hF_Error,'erf_error',saveOpts);end  
   
    % Cloud centre
    hF_Centre=showAtomCentre(erf_data,pco_xVar,ErfPopts);    
    if doSave;saveFigure(hF_Centre,'erf_position',saveOpts);end 
        
    % erf Size
    hF_size=showSize(erf_data,pco_xVar,ErfPopts);    
    if doSave;saveFigure(hF_size,'erf_size',saveOpts);end
    
    % Aspect ratio
    hF_ratio=showAspectRatio(erf_data,pco_xVar,ErfPopts);    
    if doSave;saveFigure(hF_ratio,'erf_ratio',saveOpts);end   

    % Peak erfian density
    hF_density=showDensity(erf_data,pco_xVar,ErfPopts);    
    if doSave;saveFigure(hF_density,'erf_density',saveOpts);end 
      
end