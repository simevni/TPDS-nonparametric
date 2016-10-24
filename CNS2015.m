function [betaR,rejectNull] = CNS2015(cashoffer,x,y,w,mv,c,orderPoly,nknots,alpha,beta0)
% Reference: 
% Chernozhukov, Victor, Whitney Newey, and Andres Santos (2015). "Constrained conditional moment restriction models" (No. CWP59/15). Centre for Microdata Methods and Practice, Institute for Fiscal Studies.
%-----------
%Inputs:
% cashoffer: integer or real vector (n). data. the values of x before normalization. 
% x: integer or real vector (n). data. control variable.
% y: integer or real vector (n). data. dependent variable.
% w: real vector (n). weight.
% mv: real (1). fiscal value of the rations. will be used in nonlinear
% constraint.
% c: real (1). nonlinear constraint function value.
% orderPoly: integer (1). the order of B-spline polynomials. 
% nknots: integer (1). the number of B-spline knots 
% alpha: real (1). level of significance.
% beta0: real vector (orderPoly+1+nknots). initial guess for beta.
%-----------
%Outputs:
% betaR: real vector (orderPoly+1+nknots). coefficient estimates.
% rejectNull: integer (1). 1 if the null hypothesis is rejected. 0
% otherwise.

    %% Construction of bases
    n = length(x);
    x_c = unique(x);
    nx_c = length(x_c);
		
    [P_c,knots,dP_c] = makeBsplineBasis(orderPoly,nknots,[],x_c); 
    [P,~] = makeBsplineBasis(orderPoly,nknots,knots,x); 
    
    %% A matrix for inequality constraints: Ax <= b
    A = [-dP_c; -P_c; P_c];
    b = [zeros(nx_c,1); zeros(nx_c,1); ones(nx_c,1)];
	
	ineqc_LHS = A;
	ineqc_RHS = b;

%	c = arg2(1,1); %nonlinear equality constraint 
%	mv = arg2(2,1); %0-mv is the lowest possible amount of cash offer  
%	cashoffer = arg2(:,2); %observed amount of cash offer
%	P = arg2(:,3:end); %bases of normalized observed cash offer space
	arg1 = [y w P];
	arg2 = [[[c; mv]; zeros(n-2,1)] cashoffer P];


%% Find an initial guess that satisfies the constraints

    [~,cval,otherOutputs] = nonlconfunc(beta0,arg2);   

    if (abs(cval) > 1.0e-2) %gain is too small. need to value cash more. increase the slope
        [const,fval]= fzero(@(const)nonlconfunc_aux(const,beta0(2:end),arg2),beta0(1));
        initialGuess = [const; beta0(2:end)];
    else
        initialGuess = beta0;
    end

%% Compute the test statistic I
	options = optimoptions('fmincon','Display','iter','Algorithm','sqp','Display','none');
	problem.options = options;
	problem.solver = 'fmincon';
	problem.objective = @(theta)testStatistic(theta,arg1);
	problem.x0 = initialGuess;
	problem.Aineq = ineqc_LHS; %monotonicity constraints
	problem.bineq = ineqc_RHS; 
	problem.nonlcon= @(theta)nonlconfunc(theta,arg2); %consumer surplus
	
	[betaR,testStat,exitflag,output] = fmincon(problem);	
    c1 = [ineqc_LHS*betaR ineqc_RHS];
	cdf = P*betaR;
    [~,cval,otherOutputs] = nonlconfunc(betaR,arg2); % check nonlinear constraint   
%     cval
%     otherOutputs

    
    %% Compute the critical value
    rho = y - P*betaR;
    arg2 = [[[c; mv]; zeros(n-2,1)] cashoffer P];
    cv = criticalValue (betaR,rho,P,w,ineqc_LHS,ineqc_RHS,arg2,alpha);	

    %% Check whether c is included in cv
    rejectNull = 1;
    if (testStat<=cv) 
       rejectNull = 0;
    end
    
end	%function CNS2015

%%%%%%
function [B,knots,dB] = makeBsplineBasis(p,nknots,knots,x)
% knots are created to evenly partition x 
%Input:
% p: the order of the polynomial (integer)
% m: number of knots (integer)
%knots: knots (real vector, optional)
%x: data (real vector)
%Output:
%B: the matrix of the basis
%knots: 
%dB: the matrix of the first derivative of the basis
    n = length(x); 

    if isempty(knots)==1
        %% create knots
        divSize = round(n/(nknots+1));
        mod = n - (nknots+1)*divSize;
        
        i = divSize;
        j = 1;
        while i < n - mod
            knots(j) = x(i);
            i = i + divSize;
            j = j + 1;
        end
    end

    %% Compute the basis & derivative
    B = zeros(n,p+1+nknots);
    dB = zeros(n,p+1+nknots);
    
    for i=1:p+1
        j = i-1;
        B(:,i) = x.^j; 
        dB(:,i) = j*x.^(j-1);
    end

    xLessknots = repmat(x,1,nknots)-repmat(knots,n,1);
    include = (xLessknots>=0);

    B(:,p+2:p+1+nknots) = (xLessknots.^p).*include;
    dB(:,p+2:p+1+nknots) = p*(xLessknots).^(p-1).*include;

end
	
%%%%%%	
function fval = testStatistic(beta,arg) 
% beta: the minimizer
% arg: the rest of the input arguments

	y = arg(:,1);
	w = arg(:,2); %weights
	P = arg(:,3:end);
	J = size(P,2);
	n = length(y);
	fval = norm(sum(repmat(w.*(y - P*beta),1,J).*P,1)/sqrt(n));
end %testStatistic

%%%%%%
function [ineqc,fval,otherOutputs] = nonlconfunc(theta,arg)

	c = arg(1,1); %nonlinear equality constraint 
	mv = arg(2,1); %0-mv is the lowest possible amount of cash offer  
	cashoffer = arg(:,2); %observed amount of cash offer
	
	P = arg(:,3:end); %bases of normalized observed cash offer space
	
	% compute the cdf
	cdf = P*theta;
	
	cdf = [0; cdf; 1]; 
	cashoffer = [-mv; cashoffer; cashoffer(end)*1.05];
	
	% gain from cash per hh
	gain_hh = trapz(cashoffer(cashoffer<=0),cdf(cashoffer<=0));

	% loss from cash per hh
	loss_hh = trapz(cashoffer(cashoffer>=0),1-cdf(cashoffer>=0));
	
	fval = gain_hh - loss_hh - c;
     fval = gain_hh - c;
%      fval = loss_hh - c;
	otherOutputs = gain_hh;
	
	ineqc =[]; %inequality constraint is left empty
	
end

%%%%%%
function fval = nonlconfunc_aux(const,betahat,arg)
    theta = [const; betahat];
	c = arg(1,1); %nonlinear equality constraint 
	mv = arg(2,1); %0-mv is the lowest possible amount of cash offer  
	cashoffer = arg(:,2); %observed amount of cash offer
	
	P = arg(:,3:end); %bases of normalized observed cash offer space
	
	% compute the cdf
	cdf = P*theta;
	
	cdf = [0; cdf; 1]; 
	cashoffer = [-mv; cashoffer; cashoffer(end)*1.05];
	
	% gain from cash per hh
	gain_hh = trapz(cashoffer(cashoffer<=0),cdf(cashoffer<=0));

	% loss from cash per hh
	loss_hh = trapz(cashoffer(cashoffer>=0),1-cdf(cashoffer>=0));
	
	fval = gain_hh - loss_hh - c;
     fval = gain_hh - c;
%     fval = loss_hh - c;
end
%%%%%%

function ub = criticalValue (betahat,rho,P,weight,ineqc_LHS,ineqc_RHS,arg2,alpha)
%alpha: level of significance

	n = length(rho);
	nsimul = 200;
	J = size(P,2);
	
	% Draw Wi ~ N(0,1)
	rng(29847503,'twister');
	Wi = normrnd(0,1,[n nsimul]);
	
	% Compute Un
    Un = zeros(nsimul,1);
    
    for wi=1:nsimul
        
        Gn = sum(repmat(Wi(:,wi),1,J).*...
                (repmat(weight,1,J).*P.*repmat(rho,1,J)...
                - repmat(mean(repmat(weight,1,J).*P.*repmat(rho,1,J),1),n,1)),1)/sqrt(n);

        % Compute H
        H = zeros(J);
        for i=1:n
            H = H + P(i,:)'*P(i,:);
        end
        H = H/n;

        %% Solve constrained minimization
        initialGuess = betahat;
        options = optimoptions('fmincon','Display','iter','Algorithm','sqp','Display','none');
        problem.options = options;
        problem.solver = 'fmincon';
         problem.objective = @(theta)norm(Gn'-H*theta);
%         problem.objective = @(theta)norm(Gn-mean(P.*repmat(P*theta,1,J),1));
        problem.x0 = initialGuess;
        problem.Aineq = ineqc_LHS; %monotonicity constraints
        problem.bineq = ineqc_RHS; 
        problem.nonlcon= @(theta)nonlconfunc(theta,arg2); %consumer surplus

        [~,fval] = fmincon(problem);
    	Un(wi) = fval;
        
    end %bootstrap wi
    
	%% 1-alpha percentile value
    Un = sort(Un);
    ub = Un(round(nsimul*(1-alpha)));

    
end
