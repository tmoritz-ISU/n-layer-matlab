function [Smn] = calculate(self,f,er,ur,thk)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


arguments
    self nLayerFilled
    f(:,1); 
    er(1,:); 
    ur(1,:);
    thk(1,:) {mustBeNonempty}; 
end

%% Validate Structure
[er, ur, thk] = nLayer.validateStructure(er, ur, thk, ...
    CheckStructureValues=self.checkStructureValues);

%% Pre-compute Constants
k0 = 2*pi*f/self.speedOfLight
k1 = k0.*sqrt(self.waveguidePort1er.*self.waveguidePort1ur); %k inside WG 1
kc1 = self.mode_kc0./(sqrt(self.waveguidePort1er.*self.waveguidePort1ur)) %cutoff wavenumber in WG 1
kz1 = sqrt(k1.^2-kc1.^2); %kz inside port 1
eta1 = 120*pi.*sqrt(self.waveguidePort1ur./self.waveguidePort1er); %wave impedance in WG 1

k2 = k0.*sqrt(self.waveguidePort2er.*self.waveguidePort2ur); %k inside WG 2
kc2 = self.mode_kc0./(sqrt(self.waveguidePort2er.*self.waveguidePort2ur)) %cutoff wavenumber in WG 2
kz2 = sqrt(k2.^2-kc2.^2); %kz inside WG 2
eta2 = 120*pi.*sqrt(self.waveguidePort2ur./self.waveguidePort2er); %wave impedance in WG 2

k = cell(size(thk));
kz = cell(size(thk));
eta = cell(size(thk));
Z = cell(size(thk));
Tau = cell(size(thk));
Gamma = cell(length(thk)+1, 1);

for ii=1:length(thk)
    k{ii} = k0.*sqrt(ur{ii}.*er{ii});
    kz{ii} = sqrt(k{ii}.^2-self.mode_kc0.^2);
    eta{ii} = 120*pi.*sqrt(ur{ii}./er{ii});
    Tau{ii} = exp(-1j.*kz{ii}.*thk{ii});
    if strcmp(self.modeType,"TE")
        Z{ii} = k{ii}.*eta{ii}./kz{ii};
    end
    if strcmp(self.modeType, "TM")
        Z{ii} = kz{ii}.*eta{ii}./k{ii};
    end
end

if strcmp(self.modeType,"TE")
    Z = {k1.*eta1./kz1, Z{:}, k2.*eta2./kz2}.';
end
if strcmp(self.modeType, "TM")
    Z = {kz1.*eta1./k1, Z{:}, kz2.*eta2./k2}.';
end

for ii=1:length(Z)-1
    Gamma{ii} = (Z{ii+1}-Z{ii})./(Z{ii+1}+Z{ii});
end

%% Recursively Calculate Smn
S11 = cell(size(Gamma));
S12 = cell(size(Gamma));
S21 = cell(size(Gamma));
S22 = cell(size(Gamma));

%% Initial Values
S11{1} = Gamma{1};
S12{1} = 1-Gamma{1};
S21{1} = 1+Gamma{1};
S22{1} = -Gamma{1};

for ii=1:length(thk)
    S11{ii+1} = S11{ii}+((S12{ii}.*S21{ii}.*Gamma{ii+1}.*Tau{ii}.^2)...
        ./(1-Gamma{ii+1}.*Tau{ii}.^2.*S22{ii}));
    S12{ii+1} = (Tau{ii}.*(1-Gamma{ii+1}).*S12{ii})...
        ./(1-Tau{ii}.^2.*S22{ii}.*Gamma{ii+1});
    S21{ii+1} = (S21{ii}.*Tau{ii}.*(1+Gamma{ii+1}))...
        ./(1-S22{ii}.*Tau{ii}.^2.*Gamma{ii+1});
    S22{ii+1} = -Gamma{ii+1}+((1-Gamma{ii+1}).*(1+Gamma{ii+1}).*S22{ii}.*Tau{ii}.^2)...
        ./(1-Gamma{ii+1}.*Tau{ii}.^2.*S22{ii});
end

Smn = reshape([S11{end}, S12{end}, S21{end}, S22{end}], length(f), 2, 2);

end