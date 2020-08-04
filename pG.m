% pG() - A function to plot Gaussian curve for the demo of 'Hyperbolic power transformation'
%
function [h2 xmin xmax maxfx]= pG(x,plotstyle,rotate)
  %
  mu=0;
  sigma=1;
  if nargin<3,
     rotate =0; % rotate x and y for plotting Gaussian curve on the right hand side of a figure
  end;
  %plotstyle='-';

  XMIN=min(x)-1.0; XMAX=max(x)+1.0; PNTS=20001;
  BND=max(abs(XMAX), abs(XMIN));
  x=linspace(-BND,BND,PNTS);

  f = @(x) ...
    1./(sqrt(2*pi)*sigma) .* ...
    exp( -1./(2*sigma.^2).*(x-mu).^2);

  %idxsmall=find(f(x)< 0.01);
  %x(idxsmall)=[];

  if rotate, h2=plot(f(x), x, plotstyle); 
  else, h2=plot(x, f(x), plotstyle);
  end;
  maxfx=max(f(x));
  xmin=min(x); xmax=max(x);
  %xlabel('y');

