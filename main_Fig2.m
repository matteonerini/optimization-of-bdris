% Compute received power with analytical solution.
% * SISO systems with direct link.
% * Single, group, and fully connected RISs.

clear; clc;
rng(3);
tic;

% Parameters
nMonte = 500; %5000
NIs = 1:66; %1:66;
NGs = [1,2,3,4,0]; %[1,2,3,4,0];
PT = 10; % Transmit power [W]
K = 0; % Rician factor

% Main loop
PR = nan(nMonte,length(NIs),length(NGs),2);
PR_noRIS = nan(nMonte,1);
for iMonte = 1:nMonte
    if mod(iMonte,100) == 0
        fprintf(['iMonte: ',num2str(iMonte),'\n'])
    end

    for iNI = 1:length(NIs)
        NI = NIs(iNI);
        
        % Generate channels hRT, hIT and hRI
        [GRT,GRI,GIT] = func_path_gain();
        hRT_LoS = exp(1i * 2 * pi * rand(1));
        hRT_NLoS = sqrt(1/2) * (randn(1) + 1i * randn(1));
        hRT = sqrt(GRT) * (sqrt(K/(1+K)) * hRT_LoS + sqrt(1/(1+K)) * hRT_NLoS); % Rician
        hRI_LoS = exp(1i * 2 * pi * rand(1,NI));
        hRI_NLoS = sqrt(1/2) * (randn(1,NI) + 1i * randn(1,NI));
        hRI = sqrt(GRI) * (sqrt(K/(1+K)) * hRI_LoS + sqrt(1/(1+K)) * hRI_NLoS); % Rician
        hIT_LoS = exp(1i * 2 * pi * rand(NI,1));
        hIT_NLoS = sqrt(1/2) * (randn(NI,1) + 1i * randn(NI,1));
        hIT = sqrt(GIT) * (sqrt(K/(1+K)) * hIT_LoS + sqrt(1/(1+K)) * hIT_NLoS); % Rician

        % Normalize channels hIT and hRI
        hRI_norm = hRI / norm(hRI);
        hIT_norm = hIT / norm(hIT);

        for iNG = 1:length(NGs)
            NG = NGs(iNG);

            if mod(NI,NG) == 0 || NG == 0
                % Compute Theta
                Theta = func_theta(hRI_norm,hIT_norm,NG);
                Theta = exp(1i * angle(hRT)) * Theta;
                
                % Compute received signal power
                PR(iMonte,iNI,iNG,1) = PT * abs(hRT + hRI*Theta*hIT) ^ 2;
                PR(iMonte,iNI,iNG,2) = PT * (sqrt(func_upper_bound_GC(hIT, hRI, NG)) + abs(hRT)) ^ 2;
            end
        end
    end
    PR_noRIS(iMonte) = PT * abs(hRT) ^ 2;
end

PR_av = squeeze(mean(PR));
PR_noRIS_av = mean(PR_noRIS);
toc;

%% Plot
figure('DefaultAxesFontSize',12); hold on;
LineW = 1.5;
MarkS = 8;

plot(NIs(8:8:64),PR_av(8:8:64,5,2)*1e6,'-h','linewidth',LineW,'MarkerSize',MarkS)
plot(NIs(8:8:64),PR_av(8:8:64,5,1)*1e6,'--p','linewidth',LineW,'MarkerSize',MarkS)
plot(NIs(8:8:64),PR_av(8:8:64,4,2)*1e6,'-v','linewidth',LineW,'MarkerSize',MarkS)
plot(NIs(8:8:64),PR_av(8:8:64,4,1)*1e6,'--^','linewidth',LineW,'MarkerSize',MarkS)
plot(NIs(12:6:66),PR_av(12:6:66,3,2)*1e6,'->','linewidth',LineW,'MarkerSize',MarkS)
plot(NIs(12:6:66),PR_av(12:6:66,3,1)*1e6,'--<','linewidth',LineW,'MarkerSize',MarkS)
plot(NIs(8:8:64),PR_av(8:8:64,2,2)*1e6,'-s','linewidth',LineW,'MarkerSize',MarkS)
plot(NIs(8:8:64),PR_av(8:8:64,2,1)*1e6,'--d','linewidth',LineW,'MarkerSize',MarkS)
plot(NIs(8:8:64),PR_av(8:8:64,1,2)*1e6,'-*','linewidth',LineW,'MarkerSize',MarkS)
plot(NIs(8:8:64),PR_av(8:8:64,1,1)*1e6,'--o','linewidth',LineW,'MarkerSize',MarkS)

grid on;
if K == 0
    title('Rayleigh fading')
else
    title('Rician fading - Rician factor = 3 dB')
end
xlabel('Number of RIS elements');
ylabel('Average received signal power [uW]')
legend('FC - Upper Bound','Alg. 1',...
       'GC (Group Size 4) - Upper Bound','Alg. 1',...
       'GC (Group Size 3) - Upper Bound','Alg. 1',...
       'GC (Group Size 2) - Upper Bound','Alg. 1',...
       'SC - Upper Bound','Alg. 1',... %'No RIS',...
       'location','northwest','numColumns',2);

plots=get(gca, 'Children');
legend(plots([10,8,6,4,2,9,7,5,3,1]));
ax = gca;
ax.XTick = 0:8:64;
ax.XLim = [0 64];
ax.YTick = 0:0.2:1.2;
set(gcf, 'Color', [1,1,1]);
set(gca, 'LineWidth',1);