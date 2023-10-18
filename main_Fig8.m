% Plot the complexity

clear; clc;

NI = 8:8:64;

NG = [NI;
      ones(1,length(NI)) * 4;
      ones(1,length(NI)) * 3;
      ones(1,length(NI)) * 2];

% Quasi-Newton
qn = zeros(size(NG));
for i = 1:length(NI)
    qn(:,i) = (NI(i)*(NG(:,i)+1)/2) .^ 2;
end

% Upper Bound
ub = zeros(size(NG));
for i = 1:length(NI)
    ub(:,i) = NG(:,i).^2 * NI(i);
end

%% Plot
figure('DefaultAxesFontSize',12);
LineW = 1.5;
MarkS = 8;

semilogy(NI,qn(1,:),'-h','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','FC - QN')
hold on;
grid on;
semilogy(NI,qn(2,:),'-v','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','GC (Group Size 4) - QN')
semilogy(NI,qn(3,:),'->','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','GC (Group Size 3) - QN')
semilogy(NI,qn(4,:),'-s','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','GC (Group Size 2) - QN')
set(gca,'ColorOrderIndex', 1);
semilogy(NI,ub(1,:),'--h','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Alg. 1')
semilogy(NI,ub(2,:),'--v','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Alg. 1')
semilogy(NI,ub(3,:),'-->','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Alg. 1')
semilogy(NI,ub(4,:),'--s','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Alg. 1')

xlabel('Number of RIS elements');
ylabel('Computational complexity');
legend('location','northwest','NumColumns',2);

ax = gca;
ax.XTick = 0:8:64;
ax.XLim = [0 64];
ax.YLim = [10 10^7];
set(gcf, 'Color', [1,1,1]);
set(gca, 'LineWidth',1);