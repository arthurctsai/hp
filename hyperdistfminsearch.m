% hyperdistfminsearch() - An MLE for 'Hyperbolic power transformation'
%
function [alpha, betaminus, lambdaminus, betaplus, lambdaplus,l] = hyperdistmle(x_orig, alpha, betaminus, lambdaminus, betaplus, lambdaplus)

verbos=0;

xsize = size(x_orig);

if(xsize(1) ~= 1),
    x_orig = x_orig';
end;

k=1;

% ==== i==1, for x >=0; i==2 for x<0 ====
psioveralpha=0;
for i=1:2,

  if(i==1), % for x<0, side='-'
    x = x_orig(x_orig <  0); 
    Lambdaminus(1) = lambdaminus;
    Betaminus(1) = betaminus;
    beta=betaminus;
    lambda=lambdaminus;
  else, % for x>=0, side='+'
    x = x_orig(x_orig >= 0); 
    Lambdaplus(1) = lambdaplus;
    Betaplus(1) = betaplus;
    beta=betaplus;
    lambda=lambdaplus;
  end;

  parameters0 = [beta, lambda];
  % print out the initial beta, lambda and its corresponding likelihood
  if(i==1), % for x<0, side='-'
      if verbos, disp(['- k=0 beta=' num2str(beta) ' lambda='  num2str(lambda) ' l=' num2str(likelihood(parameters0, x,alpha)) ]); end;
  else,
      if verbos, disp(['+ k=0 beta=' num2str(beta) ' lambda='  num2str(lambda) ' l=' num2str(likelihood(parameters0, x,alpha)) ]); end;
  end;


  [parameters l] = fminsearch(@(parameters) likelihood(parameters, x,alpha), parameters0);
  % to see the usage of fminsearch plz ref. http://www.mathworks.com/help/matlab/ref/fminsearch.html

  beta = parameters(1); lambda=parameters(2);
     
  betax=beta.*x;
  psioveralpha = [psioveralpha  sinh(betax).*((sech(betax)).^lambda)/beta];

  k = k + 1;

  if(i==1), % for x<0, side='-'
    if verbos, disp(['- k=' num2str(k-1) ', beta=' num2str(beta) ' lambda='  num2str(lambda) ' l=' num2str(l) ]); end;
        Lambdaminus(k)=lambda;
        Betaminus(k)=beta;
        Lminus(k)=l;
    else, % for x>=0, side='+'
    if verbos, disp(['+ k=' num2str(k-1) ', beta=' num2str(beta) ' lambda='  num2str(lambda) ' l=' num2str(l) ]); end;
        Lambdaplus(k)=lambda;
        Betaplus(k)=beta;
        Lplus(k)=l;
    end;

    if(i==1), % for x<0, side='-'
      lambdaminus = lambda;
      betaminus = beta;
    else, % for x>=0, side='+'
      lambdaplus = lambda;
      betaplus = beta;
    end;

end;

alpha = ( sum(psioveralpha.^2)./length(x_orig) ).^(-0.5);


function l = likelihood(parameters,x,alpha)
% Here FMINSEARCH is used. Within the objective function, it 
% checks to see if the parameters are within 
% the bounds. If not, return a large value for the objective function.
  beta = parameters(1); lambda=parameters(2);
  if (lambda <=1) && (beta >0) && (lambda >=-5),
      betax=beta.*x;
      psi=(alpha/parameters(1)).*sinh(betax).*((sech(betax)).^parameters(2));
      psi2=(1./parameters(1)).*sinh(betax).*((sech(betax)).^parameters(2));

      l   = -(- sum( psi.^2 ) + sum( log( alpha.*(1-lambda.*tanh(betax).^2).* sech(betax).^(lambda-1) ) ) );
      %l   = sum( psi2.^2 );
  else,
    l=1000000000;
  end;

