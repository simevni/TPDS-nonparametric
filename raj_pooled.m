    % MATLAB code to estimate the distribution with monotonicity constraint using Rajasthan stated preference data
% using the method proposed by Horowitz and Lee (2015)

% References

%% Set up
cd 'D:\Dropbox\TPDS\Choice Experiment Paper\Data\Constructed\'
% number of beneficiary households in Raj
nhh_bpl = 2.6; %million
nhh_apl = 7.4;
nhh = nhh_bpl + nhh_apl;


data = dlmread('Raj_pds_value.csv',',',1,0);

pdsvalue = data;

% %% Histogram
% figure;
% h = histogram(pdsvalue,20);
% h.FaceColor = 'w';
% h.EdgeColor = 'k';
% axis([-2500 1500 0 60]);
% ax1 = gca; % current axes
% % bar(bar_x,bar_y)
% % xlabel('Amount of cash offer (deviation from FV)');
% % ylabel('Cash take-up rates');
% % 
% % % title('Empirical Distribution of CEV in Rajasthan (All Households)');
% % saveas(gcf,'D:\Dropbox\TPDS\Choice Experiment Paper\Output\Nonparametric\Rajasthan\Raj empirical CDF diff.png');
% 
% %% Kernel Density Estimation
% [pdf,x] = ksdensity(pdsvalue);
% 
% ax1_pos = ax1.Position; % position of first axes
% ax2 = axes('Position',ax1_pos,...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','none');
% hold on
% plot(x,pdf,'Parent',ax2,'LineWidth',1.5);

%% Construct CDF from the density
[cdfx,x] = ksdensity(pdsvalue,pdsvalue,'function','cdf');
x0 = x;
cdf0 = cdfx;

x = x0;
cdfx = cdf0;

% Compute the Consumer Surplus
bin = x;
ind0 = (bin<=0);
nbin0 = sum(ind0);
nbin = length(bin);

% gain from cash per hh
gain_hh = trapz(bin(1:nbin0),cdfx(1:nbin0))

% loss from cash per hh
loss_hh = trapz(bin(nbin0+1:nbin),1-cdfx(nbin0+1:nbin))


alpha = 0.025;
nboot=100;
n = length(pdsvalue);
cdfmat = zeros(length(cdfx),nboot);
rng(85768939,'twister');

choicegain = zeros(nboot,1);
cashloss = zeros(nboot,1);
cashnetgain = zeros(nboot,1);


for i=1:nboot
    sample = sort(datasample(pdsvalue,n));
    [cdfx,x] = ksdensity(sample,sample,'function','cdf');
    duplicated = ([1; diff(x)]==0);
    small = x0(x0 < x(1));
    large = x0(x0 > x(end));
    cdfmat(:,i) = [cdf0(x0 < x(1));...
                   interp1(x(duplicated==0),cdfx(duplicated==0),x0(x0>=x(1) & x0<=x(end)));...
                   cdf0(x0>x(end))];
               
   %% 95% CI for the welfare gain/loss
    x = x0;
    cdfx = cdfmat(:,i);
   
    % Compute the Consumer Surplus
    bin = x;
    ind0 = (bin<=0);
    nbin0 = sum(ind0);
    nbin = length(bin);

    % gain from cash per hh
    gain_hh = trapz(bin(1:nbin0),cdfx(1:nbin0));

    % loss from cash per hh
    loss_hh = trapz(bin(nbin0+1:nbin),1-cdfx(nbin0+1:nbin));
    
    choicegain(i,1) = gain_hh;
    cashloss(i,1) = loss_hh;
    cashnetgain(i,1) = gain_hh - loss_hh;
   
end

choicegain = sort(choicegain);
cashloss = sort(cashloss);
cashnetgain = sort(cashnetgain);

ci95_choicegain = choicegain(round(alpha*nboot)); % 95% point-wise confidence band
ci95_choicegain(2) = choicegain(round((1-alpha)*nboot))';

ci95_cashloss = cashloss(round(alpha*nboot)); % 95% point-wise confidence band
ci95_cashloss(2) = cashloss(round((1-alpha)*nboot))';

ci95_cashnetgain = cashnetgain(round(alpha*nboot)); % 95% point-wise confidence band
ci95_cashnetgain(2) = cashnetgain(round((1-alpha)*nboot))';

sorted_cdf = sort(cdfmat,2);
ci95_cdf = sorted_cdf(:,round(alpha*nboot)); % 95% point-wise confidence band
ci95_cdf(:,2) = sorted_cdf(:,round((1-alpha)*nboot))';

x = x0;
cdfx = cdf0;

figure;
subplot(2,1,1) %add cdf in (nbw x 2) grid
plot(x,cdfx,'-k','LineWidth',1.5); %cdf
hold('on');
plot(x,ci95_cdf(:,1),':k',x,ci95_cdf(:,2),':k','LineWidth',1) % 95% confidence band
title('CDF of Valuation of in-Kind Transfers');
ylabel('cdf');
xlabel('Deviation from FV');
axis([x(1) x(end) 0 1]);

% Demand curve
demand_x = 1-cdfx; 
demand_y = x;
subplot(2,1,2) %add cdf in (nbw x 2) grid
plot(demand_x,demand_y,'k-',demand_x,zeros(size(demand_x)),'k:','LineWidth',1.5);
xlabel('Fraction of households');
ylabel({'Cash amount', '(Deviation from FV)'});
title('Demand for in-Kind Transfers over Cash');


    


