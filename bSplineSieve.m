function [betaR,P,knots] = bSplineSieve(x,y,w,orderPoly,nknots)
% Find the coefficients that minimizes the predicted and the target
% outcomes using linear sieve estimation method with monotonicity constraint.
% B-spline basis are used.
%Input
% x: integer or real vector (n). control variable in the data. x may include duplicated values.
% y: integer or real vector (n). target outcome.
% w: integer or real vector (n). weight (optional)
% orderPoly: integer (1). the order of B-spline polynomials.
% nknots: integer (1). number of knots for B-spline basis
%Output
% betaR: coefficient estimates
% P: integer or real matrix (n,orderPoly+1+nknots). regression basis
% knots: integer or real vector (nknots).


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

