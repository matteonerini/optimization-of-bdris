function [GRT,GRI,GIT] = func_path_gain_transmissive_mode()
% Generate the path gains in a single-user RIS-aided MIMO system

% Variables
locT = [0;0];  % Tx location [m]
locI = [50;2]; % RIS location [m]
locR = [52;4]; % Rx location [m]
C0 = 1e-3;     % Pathloss at D0
D0 = 1;        % Reference distance [m]
aRT = 4;       % Pathloss exponent for RT
aRI = 2.8;     % Pathloss exponent for RI
aIT = 2.0;     % Pathloss exponent for IT

% Compute the path gains
dRT = norm(locR - locT);
GRT = C0 * (dRT / D0) ^ (-aRT);

dRI = norm(locR - locI);
GRI = C0 * (dRI / D0) ^ (-aRI);

dIT = norm(locI - locT);
GIT = C0 * (dIT / D0) ^ (-aIT);

end