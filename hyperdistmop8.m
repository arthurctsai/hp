function [a, betaminus, lambdaminus, betaplus, lambdaplus, x_sorted, x_unsorted] = hyperdistmop8(x,varargin)
% 2012.08.31 Arthur at UCSD
% A newly developed hyper transformation with scalar to extend 
% the ability of the transformation.
% May 31, 2013 

                 %  name          type     range       default
defaultsetting = {'verbos'       ,'string', '',        'off' ;...
                  'visible'      ,'string', '',        'off' ;...
                  'stdx'         ,'string', '',        'on' ;...
                  'rescaling'    ,'string', '',        'off' }; 
              
g = finputcheck(varargin,defaultsetting);

varargintmp=fieldnames(g);
for i=1:length(varargintmp)
  eval([varargintmp{i} '=getfield(g,''' varargintmp{i}  ''');']);
end

% x_orig=x; % it didn't used
[m, n]=size(x);  if m<n, x=x'; end; [n, K]=size(x); 

% ===== parameter =====

quan_h=0.95; % origin 0.95
quan_l=0.75; % origin 0.75

% ===== std ======
if strcmp(stdx, 'on'),
    MM=repmat(mean(x),n,1); 
    STD=repmat(std(x),n,1);
    x = (x-MM)./STD;
end;

% ==== median removal ====
MED=repmat(median(x),n,1);
x = x - MED;

x_unsorted = x;

% ==== trim extremes ====
%x(find(x>5))=[];
%x(find(x<-5))=[];


% ==== rescale ====
if strcmp(rescaling, 'on'),
    %2012-09-01 if boundary of x is not ouside of [-2 2], we need to rescale it
    % in order to calculate Z_q and Z_p
    minboundary=min(abs(min(x)), abs(max(x)));
    ratio = 2.5./minboundary;

    [n, K]=size(x);
    ratio_tmp(1:K)=1;ratio_tmp(ratio>1)=ratio(ratio>1);ratio_tmp=repmat(ratio_tmp,n,1);
    x=ratio_tmp.*x;

    % do the same things for the unsorted x;
    [n_unsorted, K_unsorted]=size(x_unsorted);
    ratio_unsorted(1:K_unsorted)=1;ratio_unsorted(ratio>1)=ratio(ratio>1);ratio_unsorted=repmat(ratio_unsorted,n_unsorted,1);
    x_unsorted=ratio_unsorted.*x_unsorted;
end;


MED=repmat(median(x),n,1); x = x - MED;

% ==== sort ====
x=sort(x);
x_sorted=x;


% =========================================
% ==== estimate betaminus and betaplus ====
  
% negative
  side='-';%   sgn=-1;
  
  quan_hv=quantile(x,1-quan_h);
  quan_lv=quantile(x,1-quan_l);
  
  qtmpsgn=abs(quan_hv/2) < abs(quan_lv);

  qtmpxqT=(1-quan_h).*ones(1,K);  
  [~ ,qtmpxpT]=min(abs(quan_hv(ones(1,n),:)/2-x));
  qtmpxpT=qtmpxpT/n;
  xpT=quan_hv/2;
  
  qtmpxpF=(1-quan_l)*ones(1,K);  
  [~ ,qtmpxqF]=min(abs(quan_lv(ones(1,n),:)*2-x));  
  qtmpxqF=qtmpxqF/n;
  xpF=quan_lv;
  
  A=qtmpsgn.*qtmpxqT+not(qtmpsgn).*qtmpxqF; Z_q = norminv(A); % xq
  B=qtmpsgn.*qtmpxpT+not(qtmpsgn).*qtmpxpF; Z_p = norminv(B); % xp
  xp=qtmpsgn.*xpT+not(qtmpsgn).*xpF;
  
  b=zeros(1,K);
  tmp = acosh(Z_q./(2*Z_p));
  b(Z_q./(2*Z_p)>1)=tmp(Z_q./(2*Z_p)>1);

  tmp = atanh(sqrt(2*Z_p./Z_q-1));
  b(Z_q./(2*Z_p)==1)=tmp(Z_q./(2*Z_p)==1);
  b(Z_q./(2*Z_p)<1)=tmp(Z_q./(2*Z_p)<1);


if strcmp(verbos, 'on'),
  if Z_q./(2*Z_p)>1,
    disp(['for the ' side '  side Z_q./(2*Z_p)=' num2str(Z_q./(2*Z_p)) '>1 (for acosh) =>b' side '=' num2str(b)]);
  else
    disp(['for the ' side ' side 2*Z_p/Z_q-1: -1<' num2str(2*Z_p/Z_q-1) '<1 (for atanh) => b' side '=' num2str(b)]);
  end;  
end

betaminus = b./abs(xp); 

% positive

  side = '+';%   sgn = 1;
  
  quan_hv=quantile(x,quan_h);
  quan_lv=quantile(x,quan_l);
  
  qtmpsgn=abs(quan_hv/2) < abs(quan_lv);

  qtmpxqT=(quan_h).*ones(1,K);  
  [~ ,qtmpxpT]=min(abs(quan_hv(ones(1,n),:)/2-x));
  qtmpxpT=qtmpxpT/n;
  xpT=quan_hv/2;
  
  qtmpxpF=(quan_l).*ones(1,K);  
  [~ ,qtmpxqF]=min(abs(quan_lv(ones(1,n),:)*2-x));  
  qtmpxqF=qtmpxqF/n;
  xpF=quan_lv;
  
  A=qtmpsgn.*qtmpxqT+not(qtmpsgn).*qtmpxqF; Z_q = norminv(A); % xq
  B=qtmpsgn.*qtmpxpT+not(qtmpsgn).*qtmpxpF; Z_p = norminv(B); % xp
  xp=qtmpsgn.*xpT+not(qtmpsgn).*xpF;


  b=zeros(1,K);
  tmp = acosh(Z_q./(2*Z_p));
  b(Z_q./(2*Z_p)>1)=tmp(Z_q./(2*Z_p)>1);

  tmp = atanh(sqrt(2*Z_p./Z_q-1));
  b(Z_q./(2*Z_p)==1)=tmp(Z_q./(2*Z_p)==1);
  b(Z_q./(2*Z_p)<1)=tmp(Z_q./(2*Z_p)<1);


if strcmp(verbos, 'on'),
  if Z_q./(2*Z_p)>1,
      disp(['for the ' side '  side Z_q./(2*Z_p)=' num2str(Z_q./(2*Z_p)) '>1 (for acosh) =>b' side '=' num2str(b)]); 
  else
      disp(['for the ' side ' side 2*Z_p/Z_q-1: -1<' num2str(2*Z_p/Z_q-1) '<1 (for atanh) => b' side '=' num2str(b)]);
  end;
end

betaplus = b./abs(xp); 

%return;
%
% ==== estimate lambda ====

C=0.15;D=0.4;
Z=norminv(C)/norminv(D); % Z_s./Z_t

qx=quantile(x,[C D 1-C 1-D]); if size(qx,1)==1, qx=qx'; end

mylambda  = ( log(Z)+log(sinh(qx([2 4],:).*[betaminus;betaplus]))-log(sinh(qx([1 3],:).*[betaminus;betaplus]))) ...
   ./(log(sech(qx([1 3],:).*[betaminus;betaplus]))-log(sech(qx([2 4],:).*[betaminus;betaplus])) );

mylambda(mylambda>1)=1; mylambda(mylambda<-1)=-1;

lambdaminus  = mylambda(1,:);
lambdaplus  = mylambda(2,:);


if strcmp(verbos, 'on'),
   disp(['lambdaminus=' num2str(b)]); 
   disp(['lambdaplus=' num2str(b)]);
end



zDMinus=norminv(D); % Z_t
xDMinus=qx(2,:);
zDPlus=-zDMinus;    % -Z_t
xDPlus=qx(4,:);

% calculate alpha by mehtod of percentiles
a=(betaplus.*zDPlus-betaminus.*zDMinus)./(sinh(betaplus.*xDPlus).*sech(betaplus.*xDPlus).^lambdaplus - sinh(betaminus.*xDMinus).*sech(betaminus.*xDMinus).^lambdaminus);

if strcmp(visible, 'on'),
    figure; 
    %set(gcf,'Position',get(0,'ScreenSize'));
    fig_m=round(sqrt(K));
    fig_n=ceil(sqrt(K));

    for i=1:K
        subplot(fig_m,fig_n,i);
        myhist_simple(x(:,i),60); hold on;
%         [h2 xmin xmax maxfx]=px(x(:,i), a(i), betaminus(i), lambdaminus(i), betaplus(i), lambdaplus(i), 'r-');
%         axis([-5, 5,0,1.1*maxfx]);
        px(x(:,i), a(i), betaminus(i), lambdaminus(i), betaplus(i), lambdaplus(i), 'r-');
        title(['s_' num2str(i)]);
    end
end;

