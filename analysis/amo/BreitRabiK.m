function out = BreitRabiK(B,F,mF)
% out = BreitRabiK(B, F, mF)
%
% Calculates hyperfine energy for K state |F, mF> at magnetic field B (in G) 
% using the Breit-Rabi formula (source: Ramsey, "Molecular Beams", p.86; 
% here using notation from T.G. Tiecke, "Proerties of Potassium"). Output 
% in units of energy [J]. Transcribed from established Mathematica code.
% S. Trotzky, May 2014

    h = 6.6260755e-34;
    muB = 9.27400915e-24; % (J/T)
    ahfs = -h*285730800; % is ground state splitting divided by 9/2
    gS = 2.0023193043622; % electron-g (equal to "gJ" since S = J = 1/2)
    gI = 0.000176490; % 40K (Tiecke)
        
    xK = @(B) ((gS - gI)*muB*B / (9/2*ahfs));
    ebrK = @(B,F,mF) ( -ahfs/4 + gI*muB*B*1e-4.*mF + (-1).^(F-1/2)*9/2*ahfs/2*sqrt( 1 + 4*mF/9.*xK(B*1e-4) + xK(B*1e-4).^2 ) );
    
    out = ebrK(B,F,mF);
    
end

% Another function definition which avoids dealing with small numbers
function out = BreitRabiK2(B,F,mF)

% Hyperfine constant in MHz
ahf = 285.7308;

% Bohr magneton divided by (MHz/Gauss)
muB = 1.39962449; 

% Coupling constant
gS = 2.0023193043622; % electron-g (equal to "gJ" since S = J = 1/2)
gI = 0.000176490; % 40K (Tiecke)
    
% Total nuclear spin
I = 4;

% x parameter is a function of magnetic field
xB = @(B) (muB/ahf)*((gS-gI)/(I+1/2))*B;

% Nuclear zeeman coupling parameter
c1 = @(B) gI*(muB/ahf)*B;

ebrK = @(B,F,mF) ahf*(-0.25 + c1(B)*mF*B + ...
    (-1).^(F-1/2)*((I+1/2)/2)*...
    sqrt(1 + (4*mF/(2*I+1)).*xB(B) + xB(B).^2 ));

out = ebrK(B,F,mF);

end