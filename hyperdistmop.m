% hyperdistmop() - The method of percentiles for 'Hyperbolic power transformation'
%
function [a betaminus, lambdaminus, betaplus, lambdaplus, x] = hyperdistmop(x)

verbos=0; alphaByMLE='off';

if strcmp(alphaByMLE,'on'), % note that if you calculate alpha by MLE, outlier have to be removed to stabilize the variance
  % 2012.11.6 outlier removal, note that now works only for one dimentional data
  [index] = outlier(x, 0.05, length(x)*0.03);
  x(index)=[]; % maybe you can just put randomly its sample instead of get rid of them
end;


[m, n]=size(x); if m<n, x=x'; end; [n, K]=size(x); 

MM=repmat(mean(x),n,1); 
STD=repmat(std(x),n,1);
x = (x-MM)./STD;
MED=repmat(median(x),n,1);
x = x - MED;

%2012-09-01 if boundary of x not ouside of [-2 2], we need to rescale it
% in order to calculate z_A and z_B
minboundary=min(abs(min(x)), abs(max(x)));
ratio = 2.5./minboundary;
if ratio>1,
  x=ratio*x;
end;

MED=repmat(median(x),n,1); x = x - MED;

x_unsorted = x;
x=sort(x);
%n=length(x);

%keyboard;
% ==== estimate betaminus and betaplus ====
for k=1:K,
for i=1:2,
  % make side=='-' or '+'
  sgn=i-1; if ~sgn, sgn=-1; end; side=num2str(sgn); side=side(1); if strcmp(side,'1'), side='+'; end;

  [b idx]=min(abs(x(:,k)-sgn*2));
  A=idx./n; z_A = norminv(A);
  %C
  [b idx]=min(abs(x(:,k)-sgn*1));
  B=idx./n; z_B = norminv(B);
  %D
  tmp=2*z_B./z_A-1;
  if (tmp>-1) && (tmp<1),
    b=atanh(sqrt(tmp)); %keyboard;
    mylambda=1;
    if verbos, disp(['for the ' side ' side 2*z_B/z_A-1: -1<' num2str(2*z_B/z_A-1) '<1 (for atanh) => b' side '=' num2str(b)]); end;
  end;

  tmp=z_A./(2*z_B);
  if tmp > 1,
    b=acosh(tmp);
    mylambda=0;
    if verbos, disp(['for the ' side '  side z_A./(2*z_B)=' num2str(z_A./(2*z_B)) '>1 (for acosh) =>b' side '=' num2str(b)]); end;
  end;
  
  if ~exist('mylambda', 'var'), keyboard; end;

  %testing fix b to .5
  %b=.5;
  if sgn<0, 
      eval(['betaminus(' num2str(k) ')=' num2str(b) ';']); eval(['lambdaminus(' num2str(k) ')=' num2str(mylambda) ';']); 
  else
      eval(['betaplus(' num2str(k) ')=' num2str(b) ';']); eval(['lambdaplus(' num2str(k) ')=' num2str(mylambda) ';']);
  end;
end;
end;

%return;
%
% ==== estimate lambdaminus ====
C = 0.15; z_C = norminv(C); % note. 2012-09-25 C=0.15, 0.45 good for data2
D  = 0.4; z_D = norminv(D); % note. 2012-09-25 C=0.15, 0.25 good for data2
gamma = 0.35; z_gamma = norminv(gamma);

qx = quantile(x,[C D gamma]);
b=betaminus;
mylambda = ( log(z_C./z_D)+log(sinh(qx(2)*b))-log(sinh(qx(1)*b)) )./(log(sech(qx(1)*b))-log(sech(qx(2)*b)) );
mylambda(mylambda>1)=1; mylambda(mylambda<-1)=-1;
lambdaminus  = mylambda;


zDMinus=z_D;
xDMinus=qx(2);
a1 = z_gamma*b./(sinh(qx(3).*b).*(sech(qx(3)*b)).^mylambda);


% ==== estimate lambdaplus ====
C = 1-C; z_C = norminv(C);
D  = 1-D; z_D = norminv(D);
gamma = 1-0.35; z_gamma = norminv(gamma);

qx = quantile(x,[C D gamma]);
b=betaplus;
mylambda  = ( log(z_C./z_D)+log(sinh(qx(2).*b))-log(sinh(qx(1).*b)) )./(log(sech(qx(1).*b))-log(sech(qx(2).*b)) );
mylambda(mylambda>1)=1; mylambda(mylambda<-1)=-1;
lambdaplus  = mylambda;


zDPlus=z_D;
xDPlus=qx(2);
a2= z_gamma.*b./(sinh(qx(3)*b).*(sech(qx(3).*b)).^mylambda);

% calculate alpha by mehtod of percentiles
a=(betaplus.*zDPlus-betaminus.*zDMinus)./(sinh(betaplus.*xDPlus).*sech(betaplus.*xDPlus).^lambdaplus - sinh(betaminus.*xDMinus).*sech(betaminus.*xDMinus).^lambdaminus);
%keyboard

if strcmp(alphaByMLE,'on'),
  % == calculating alpha by MLE=======================
  x_orig = x;
  xsize = size(x_orig);
  if(xsize(1) ~= 1), x_orig = x_orig'; end;

    % ==== i==1, for x >=0; i==2 for x<0 ====
    psioveralpha=0;
    for i=1:2,

      if(i==1), % for x<0, side='-'
        x = x_orig(x_orig <  0); 
        beta=betaminus;
        lambda=lambdaminus;
    else, % for x>=0, side='+'
      x = x_orig(x_orig >= 0); 
      beta=betaplus;
      lambda=lambdaplus;
    end;

    betax=beta.*x;
    psioveralpha = [psioveralpha  sinh(betax).*((sech(betax)).^lambda)/beta];

  end;

  %disp('hyperdistmop.m:');
  %size(psioveralpha)

  alpha = ( sum(psioveralpha.^2)./length(x_orig) ).^(-0.5);
  %keyboard;

  a=alpha;

  x=x_orig;
end;


x=x_unsorted;
