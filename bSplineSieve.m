function [betaR,P,knots] = bSplineSieve(x,y,w,orderPoly,nknots)
%Input:
% x: x values. x may include duplicated values.
% x_c: unique values of x
% y: target outcome
% w: weight
% c: nonlinear constraint const.
% orderPoly:
% nknots: number of knots for B-spline basis
% knots_in: knots. real vector with length=nknots. optional.
% inference: 0: no test, 1: compute critical value with level of
% significance = alpha
%Output:
% betaR: coefficient estimates
% P: regression basis
% testResult: 0 if the computed test statistic < critical value, 1 otherwise.

    %% Construction of bases
    x_c = unique(x);
    nx_c = length(x_c);
		
    [P_c,knots,dP_c] = makeBsplineBasis(orderPoly,nknots,[],x_c); 
    [P,~] = makeBsplineBasis(orderPoly,nknots,knots,x); 
    	  
    %% A matrix for inequality constraints: Ax <= b
    A = [-dP_c; -P_c; P_c];
    b = [zeros(nx_c,1); zeros(nx_c,1); ones(nx_c,1)];
	
	%% Solve the constrained minimization
    
    % Find an initial guess that satisfies the constraints
    betahat = inv(P'*diag(w)*P)*P'*diag(w)*y;  
    betaR = quadprog(P'*diag(w)*P,-P'*diag(w)*y,A,b,[],[],[],[],betahat);
               
end



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
        if ( j == 0 )
            dB(:,i) = 0;
        else
            dB(:,i) = j*x.^(j-1);
        end
    end

    xLessknots = repmat(x,1,nknots)-repmat(knots,n,1);
    include = (xLessknots>=0);

    B(:,p+2:p+1+nknots) = (xLessknots.^p).*include;
    dB(:,p+2:p+1+nknots) = p*(xLessknots).^(p-1).*include;

end

