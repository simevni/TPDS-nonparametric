% using the method proposed by Horowitz and Lee (2015)

% References

%% Set up
cd 'D:\Dropbox\TPDS\Choice Experiment Paper\Data\Constructed\'

data_bpl = dlmread('Bihar_bpl_6.csv',',',1,0); %cash offer - MV
data_aay = dlmread('Bihar_aay_6.csv',',',1,0); %cash offer - MV

alpha = 0.05;
orderPoly = 3;

mv_bpl = 527+88; 
mv_aay = 861+88;



%%%%%%%%%%%%%
%% BPL 
%%%%%%%%%%%%%
encash = data_bpl(:,2);
cashoffer = data_bpl(:,1);
weight = data_bpl(:,3);
mv = mv_bpl;

y = encash;
n = length(y); %number of observations
x = (cashoffer-cashoffer(1))/(1+cashoffer(n)-cashoffer(1)); %bound x between 0 and 1

%% Find the optimal number of knots using cross-validation
% maxnknots = round(4*length(unique(x))^0.15);
% c =[];
% CV = zeros(maxnknots,1);
% 
% % Polynomial basis
% B = zeros(n,orderPoly);
% B(:,1) = 1;
% for i=1:orderPoly
%     B(:,i+1) = x.^i; 
% end
% 
% for nknots = 1:maxnknots
%     predictErrorsq = zeros(n,1);
%     for i=1:n
%         x_i = [x(1:i-1);x(i+1:end)]; %leave-one-out sample
%         y_i = [y(1:i-1);y(i+1:end)];
%         w_i = [weight(1:i-1);weight(i+1:end)]; 
%         [beta,~,knots] = bSplineSieve(x_i,y_i,w_i,orderPoly,nknots);
% 
%         xLessknots = repmat(x(i),1,nknots)-knots;
%         include = (xLessknots>0);
%         B_i =[B(i,:) (xLessknots.^orderPoly).*include];
%         predictErrorsq(i) = (y(i) - B_i*beta).^2;
%     end
%     CV(nknots) = mean(predictErrorsq);
% end
% 
% [~,nknots] = min(CV);

nknots = 6; 

%% Estimate the CDF using the optimal number of knots
[beta,P,~] = bSplineSieve(x,y,weight,orderPoly,nknots);
cdf = P*beta;

%% Compute the CS using the estimated CDF
cdf0 = [0; cdf; 1]; 
cashoffer0 = [-mv; cashoffer; cashoffer(end)*1.05];

% gain from cash per hh
gain_hh = trapz(cashoffer0(cashoffer0<=0),cdf0(cashoffer0<=0));

% loss from cash per hh
loss_hh = trapz(cashoffer0(cashoffer0>=0),1-cdf0(cashoffer0>=0));

c0 = gain_hh - loss_hh;

% Plot the cdf
% figure;
% subplot(3,1,1) %add cdf in (nbw x 3) grid
% plot(cashoffer,cdf,'-k','LineWidth',1.5); %cdf
% hold('on');

%% Compute the 95% CI of net welfare gain under the cash only system and plot
% Get curves that corresponds to the 95% CI of net welfare gain of cash
% only system

% % 95% CI lb
% c = loss_hh*1.145; %BPL lower bound = c0*57.55, BPL upper bound = c0*-80.32
% [betahat,testResult] = CNS2015(cashoffer,x,y,weight,mv,c,orderPoly,nknots,alpha,beta);
% ci95_cdf = P*betahat;
% % 
% % 95% CI ub
% c = c0*9.825; %BPL lower bound = c0*57.55, BPL upper bound = c0*-80.32
% [betahat,testResult] = CNS2015(cashoffer,x,y,weight,mv,c,orderPoly,nknots,alpha,beta);
% ci95_cdf(:,2) = P*betahat;
% 
% plot(cashoffer,ci95_cdf(:,1),':k',cashoffer,ci95_cdf(:,2),':k','LineWidth',1) % 95% confidence band
% title('CDF of Valuation of in-Kind Transfers (BPL)');
% ylabel('cdf');
% xlabel('Deviation from FV');
% axis([cashoffer(1) cashoffer(n) 0 1]);
%     
% % market value of the rationed goods by ration type
% demand_x_bpl = 1-[0; cdf; 1]; 
% demand_y_bpl = [-mv; cashoffer; cashoffer(n)*1.05]; %take up is zero for zero cash amount, take up is 1 for the amouht that is slightly larger than the max observed amount


%%%%%%%%%%%%%
%% AAY 
%%%%%%%%%%%%%
encash = data_aay(:,2);
cashoffer = data_aay(:,1);
weight = data_aay(:,3);
mv = mv_aay;

y = encash;
n = length(y); %number of observations
x = (cashoffer-cashoffer(1))/(1+cashoffer(n)-cashoffer(1)); %bound x between 0 and 1

%% Find the optimal number of knots using cross-validation
% maxnknots = round(4*length(unique(x))^0.15);
% c =[];
% CV = zeros(maxnknots,1);
% 
% % Polynomial basis
% B = zeros(n,orderPoly);
% B(:,1) = 1;
% for i=1:orderPoly
%     B(:,i+1) = x.^i; 
% end
% 
% for nknots = 1:maxnknots
%     predictErrorsq = zeros(n,1);
%     for i=1:n
%         x_i = [x(1:i-1);x(i+1:end)]; %leave-one-out sample
%         y_i = [y(1:i-1);y(i+1:end)];
%         w_i = [weight(1:i-1);weight(i+1:end)]; 
%         [beta,~,knots] = bSplineSieve(x_i,y_i,w_i,orderPoly,nknots);
% 
%         xLessknots = repmat(x(i),1,nknots)-knots;
%         include = (xLessknots>0);
%         B_i =[B(i,:) (xLessknots.^orderPoly).*include];
%         predictErrorsq(i) = (y(i) - B_i*beta).^2;
%     end
%     CV(nknots) = mean(predictErrorsq);
% end
% 
% [~,nknots] = min(CV);

nknots = 6;
%% Estimate the CDF using the optimal number of knots
[beta,P,~] = bSplineSieve(x,y,weight,orderPoly,nknots);
cdf = P*beta;

%% Compute the CS using the estimated CDF
cdf0 = [0; cdf; 1]; 
cashoffer0 = [-mv; cashoffer; cashoffer(end)*1.05];

% gain from cash per hh
gain_hh = trapz(cashoffer0(cashoffer0<=0),cdf0(cashoffer0<=0));

% loss from cash per hh
loss_hh = trapz(cashoffer0(cashoffer0>=0),1-cdf0(cashoffer0>=0));

c0 = gain_hh - loss_hh;

% Plot the cdf
% subplot(3,1,2) %add cdf in (nbw x 3) grid
% plot(cashoffer,cdf,'-k','LineWidth',1.5); %cdf
% hold('on');

%% Compute the 95% CI of net welfare gain under the cash only system and plot
% Get curves that corresponds to the 95% CI of net welfare gain of cash
% only system


% 95% CI lb
c = gain_hh*0.697; %BPL lower bound = c0*57.55, BPL upper bound = c0*-80.32
[betahat,testResult] = CNS2015(cashoffer,x,y,weight,mv,c,orderPoly,nknots,alpha,beta);
ci95_cdf = P*betahat;

% 95% CI ub
c = c0*2.175; %BPL lower bound = c0*57.55, BPL upper bound = c0*-80.32
[betahat,testResult] = CNS2015(cashoffer,x,y,weight,mv,c,orderPoly,nknots,alpha,beta);
ci95_cdf(:,2) = P*betahat;

plot(cashoffer,ci95_cdf(:,1),':k',cashoffer,ci95_cdf(:,2),':k','LineWidth',1) % 95% confidence band
title('CDF of Valuation of in-Kind Transfers (AAY)');
ylabel('cdf');
xlabel('Deviation from FV');
axis([cashoffer(1) cashoffer(n) 0 1]);
    
% market value of the rationed goods by ration type
demand_x_aay = 1-[0; cdf; 1]; 
demand_y_aay = [-mv; cashoffer; cashoffer(n)*1.05]; %take up is zero for zero cash amount, take up is 1 for the amouht that is slightly larger than the max observed amount



%% Demand curve

nbpl = 6.523;
naay = 2.501;

demand_y = linspace(-mv_aay,max(demand_y_bpl(end),demand_y_aay(end)),500); %create the common support for cash offer amount
demand_x_bpl = sort(unique(demand_x_bpl),'descend');
demand_x_aay = sort(unique(demand_x_aay),'descend');
demand_y_bpl = unique(demand_y_bpl);
demand_y_aay = unique(demand_y_aay);

demand_x_bpl_interpolated = interp1(demand_y_bpl,demand_x_bpl,demand_y);
demand_x_aay_interpolated = interp1(demand_y_aay,demand_x_aay,demand_y); 

demand_x_bpl_interpolated(demand_y<demand_y_bpl(1)) = 1;
demand_x_bpl_interpolated(demand_y>demand_y_bpl(end)) = 0;
demand_x_aay_interpolated(demand_y<demand_y_aay(1)) = 1;
demand_x_aay_interpolated(demand_y>demand_y_aay(end)) = 0;

demand_x = demand_x_bpl_interpolated * nbpl/(nbpl+naay) + demand_x_aay_interpolated * naay/(nbpl+naay);

subplot(3,1,3) %add demand in (nbw x 3) grid
plot(demand_x,demand_y,'k-',demand_x,zeros(size(demand_x)),'k:','LineWidth',1.5);
xlabel('Fraction of households');
ylabel({'Cash amount', '(Deviation from FV)'});
axis([demand_x(end) demand_x(1) demand_y(1) demand_y(end)]);
title('Demand for in-Kind Transfers over Cash Transfers');
    
saveas(gcf,'D:\Dropbox\TPDS\Choice Experiment Paper\Output\Nonparametric\Bihar\CDF_demand.png');

