% Compute received power with analytical solution.
% * SU-MIMO systems without direct link
% * Fully connected RISs.

clear; clc;
rng(3);
tic;

% Parameters
nMonte = 500; %5000
NIs = 8:8:64; %8:8:64;
NTNR = [[2,4,2,4,4];
        [1,1,2,2,4]];
PT = 10; % Transmit power [W]
K = 0; % Rician factor

% Main loop
PR = nan(nMonte,length(NIs),length(NTNR),2);
for iMonte = 1:nMonte
    if mod(iMonte,100) == 0
        fprintf(['iMonte: ',num2str(iMonte),'\n'])
    end

    for iNI = 1:length(NIs)
        NI = NIs(iNI);
        
        for iNTNR = 1:length(NTNR)
            NT = NTNR(1,iNTNR);
            NR = NTNR(2,iNTNR);

            % Generate channels hRT, hIT and hRI
            [GRT,GRI,GIT] = func_path_gain();
            HRT = zeros(NR,NT); % Zero
            HRI_LoS = exp(1i * 2 * pi * rand(NR,NI));
            HRI_NLoS = sqrt(1/2) * (randn(NR,NI) + 1i * randn(NR,NI));
            HRI = sqrt(GRI) * (sqrt(K/(1+K)) * HRI_LoS + sqrt(1/(1+K)) * HRI_NLoS); % Rician
            HIT_LoS = exp(1i * 2 * pi * rand(NI,NT));
            HIT_NLoS = sqrt(1/2) * (randn(NI,NT) + 1i * randn(NI,NT));
            HIT = sqrt(GIT) * (sqrt(K/(1+K)) * HIT_LoS + sqrt(1/(1+K)) * HIT_NLoS); % Rician
    
            % Normalize channels hIT and hRI
            [URI,~,~] = svd(HRI');
            uRI = URI(:,1);
            [UIT,~,~] = svd(HIT);
            uIT = UIT(:,1);

            % Compute Theta
            Theta = func_theta(uRI',uIT,0);
            
            % Compute received signal power
            PR(iMonte,iNI,iNTNR,1) = PT * norm(HRT + HRI*Theta*HIT) ^ 2;
            PR(iMonte,iNI,iNTNR,2) = PT * norm(HRI) ^ 2 * norm(HIT) ^ 2;
        end
    end
end

PR_av = squeeze(mean(PR));
toc;

%% Plot
figure('DefaultAxesFontSize',12); hold on;
LineW = 1.5;
MarkS = 8;

plot(NIs,PR_av(:,5,2)*1e6,'-h','linewidth',LineW,'MarkerSize',MarkS)
plot(NIs,PR_av(:,5,1)*1e6,'--p','linewidth',LineW,'MarkerSize',MarkS)
plot(NIs,PR_av(:,4,2)*1e6,'-v','linewidth',LineW,'MarkerSize',MarkS)
plot(NIs,PR_av(:,4,1)*1e6,'--^','linewidth',LineW,'MarkerSize',MarkS)
plot(NIs,PR_av(:,3,2)*1e6,'->','linewidth',LineW,'MarkerSize',MarkS)
plot(NIs,PR_av(:,3,1)*1e6,'--<','linewidth',LineW,'MarkerSize',MarkS)
plot(NIs,PR_av(:,2,2)*1e6,'-s','linewidth',LineW,'MarkerSize',MarkS)
plot(NIs,PR_av(:,2,1)*1e6,'--d','linewidth',LineW,'MarkerSize',MarkS)
plot(NIs,PR_av(:,1,2)*1e6,'-*','linewidth',LineW,'MarkerSize',MarkS)
plot(NIs,PR_av(:,1,1)*1e6,'--o','linewidth',LineW,'MarkerSize',MarkS)

grid on;
if K == 0
    title('Rayleigh fading')
else
    title('Rician fading - Rician factor = 3 dB')
end
xlabel('Number of RIS elements');
ylabel('Average received signal power [uW]')
legend('4 x 4 - Upper Bound','Alg. 1',...
       '2 x 4 - Upper Bound','Alg. 1',...
       '2 x 2 - Upper Bound','Alg. 1',...
       '1 x 4 - Upper Bound','Alg. 1',...
       '1 x 2 - Upper Bound','Alg. 1',...
       'location','northwest','numColumns',2);
plots=get(gca, 'Children');
legend(plots([10,8,6,4,2,9,7,5,3,1]));
ax = gca;
ax.XTick = 0:8:64;
ax.XLim = [0 64];
ax.YTick = 0:0.2:1.6;
ax.YLim = [0 1.6];
set(gcf, 'Color', [1,1,1]);
set(gca, 'LineWidth',1);