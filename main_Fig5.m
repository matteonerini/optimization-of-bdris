% Compute received power with analytical solution.
% * SU-MIMO systems with direct link.
% * Single, group, and fully connected RISs.

clear; clc;
rng(3);
tic;

% Parameters
nMonte = 500; %5000;
NIs = 8:8:64; %8:8:64;
NGs = [1,2,0]; %[1,2,0];
NTNR = [[2,2,4];
        [1,2,4]];
PT = 10; % Transmit power [W]
K = 0; % Rician factor

% Main loop
PR = nan(nMonte,length(NIs),length(NGs),length(NTNR),3);
parfor iMonte = 1:nMonte
    if mod(iMonte,100) == 0
        fprintf(['iMonte: ',num2str(iMonte),'\n'])
    end

    PR_iter = nan(length(NIs),length(NGs),length(NTNR),3);

    for iNI = 1:length(NIs)
        NI = NIs(iNI);

        for iNG = 1:length(NGs)
            NG = NGs(iNG);
        
            for iNTNR = 1:length(NTNR)
                NT = NTNR(1,iNTNR);
                NR = NTNR(2,iNTNR);
    
                % Generate channels hRT, hIT and hRI
                [GRT,GRI,GIT] = func_path_gain();
                HRT_LoS = exp(1i * 2 * pi * rand(NR,NT));
                HRT_NLoS = sqrt(1/2) * (randn(NR,NT) + 1i * randn(NR,NT));
                HRT = sqrt(GRT) * (sqrt(K/(1+K)) * HRT_LoS + sqrt(1/(1+K)) * HRT_NLoS); % Rician
                HRI_LoS = exp(1i * 2 * pi * rand(NR,NI));
                HRI_NLoS = sqrt(1/2) * (randn(NR,NI) + 1i * randn(NR,NI));
                HRI = sqrt(GRI) * (sqrt(K/(1+K)) * HRI_LoS + sqrt(1/(1+K)) * HRI_NLoS); % Rician
                HIT_LoS = exp(1i * 2 * pi * rand(NI,NT));
                HIT_NLoS = sqrt(1/2) * (randn(NI,NT) + 1i * randn(NI,NT));
                HIT = sqrt(GIT) * (sqrt(K/(1+K)) * HIT_LoS + sqrt(1/(1+K)) * HIT_NLoS); % Rician
                
                % Initialize g and w
                g = randn(1,NR) + 1i * randn(1,NR);
                g = g / norm(g);
                w = randn(NT,1) + 1i * randn(NT,1);
                w = w / norm(w);

                % Optimize Theta
                P_tmp = zeros(1,1e3);
                for iIter = 1:1e3
                    hRIeff = g * HRI;
                    hITeff = HIT * w;
                    hRTeff = g * HRT * w;

                    hRIeff_norm = hRIeff / norm(hRIeff);
                    hITeff_norm = hITeff / norm(hITeff);

                    T_tmp = func_theta(hRIeff_norm,hITeff_norm,NG);
                    T_tmp = exp(1i * angle(hRTeff)) * T_tmp;

                    [U,S,V] = svd(HRT + HRI*T_tmp*HIT);
                    g = U(:,1)';
                    w = V(:,1);
                    P_tmp(iIter+1) = (S(1,1))^2; % same as abs(g*hRI*T_tmp*hIT*w)^2;
                    % Stopping condition
                    if (P_tmp(iIter+1) - P_tmp(iIter))/P_tmp(iIter) < 1e-4
                        break;
                    end
                end
                Theta = T_tmp;

                % Compute received signal power
                PR_iter(iNI,iNG,iNTNR,1) = PT * norm(HRT + HRI*Theta*HIT) ^ 2;

            end
        end
    end
    PR(iMonte,:,:,:,:) = PR_iter;
end

PR_av = squeeze(mean(PR));
toc;

%% Plot
figure('DefaultAxesFontSize',12); hold on;
LineW = 1.5;
MarkS = 8;

plot(NIs,PR_av(:,1,3,1)*1e6,':p','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','4 x 4 - SC')
plot(NIs,PR_av(:,1,2,1)*1e6,':<','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','2 x 2 - SC')
plot(NIs,PR_av(:,1,1,1)*1e6,':o','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','1 x 2 - SC')
set(gca,'ColorOrderIndex',1)
plot(NIs,PR_av(:,2,3,1)*1e6,'--p','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','GC (Group Size 2)')
plot(NIs,PR_av(:,2,2,1)*1e6,'--<','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','GC (Group Size 2)')
plot(NIs,PR_av(:,2,1,1)*1e6,'--o','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','GC (Group Size 2)')
set(gca,'ColorOrderIndex',1)
plot(NIs,PR_av(:,3,3,1)*1e6,'-p','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','FC')
plot(NIs,PR_av(:,3,2,1)*1e6,'-<','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','FC')
plot(NIs,PR_av(:,3,1,1)*1e6,'-o','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','FC')

grid on;
if K == 0
    title('Rayleigh fading')
else
    title('Rician fading - Rician factor = 3 dB')
end
xlabel('Number of RIS elements');
ylabel('Average received signal power [uW]')
legend('location','northwest','numColumns',3);
ax = gca;
ax.XTick = 0:8:64;
ax.XLim = [0 64];
ax.YTick = 0:0.2:2;
ax.YLim = [0 2];
set(gcf, 'Color', [1,1,1]);
set(gca, 'LineWidth',1);