function [Rb,K]=getAtomicProperties
%GETATOMICPROPERTIES Summary of this function goes here
%   Detailed explanation goes here

amu=1.66053907E-27;

Rb=struct;
Rb.Atom='87Rb';
Rb.Gamma=2*pi*6E6;
Rb.Lambda=780E-9;
Rb.Mass=86.909180527*amu;
Rb.CrossSection=3/(2*pi)*Rb.Lambda.^2; % Simple two level cycling

K=struct;
K.Atom='40K';
K.Gamma=2*pi*6E6;
K.Lambda=767E-9;
K.Mass=39.96399848*amu;
K.CrossSection=3/(2*pi)*K.Lambda.^2; % Simple two level cycling

end

