function dx = model_LacIIPTGparams(t,xvals,params)


% Defining model variables:
dx = zeros(3,1);
I  = xvals(1); % internalized IPTG
L  = xvals(2); % LacI
C  = xvals(3); % LacI.IPTG Complex


% ------------------------------------------------------------------------
% [1] MODEL - PARAMETERS
% ------------------------------------------------------------------------
p = params;

% Fwd and rev rates of IPTG sequestration of LacI:
  kf    = p(6);  %
  kr    = p(7);  %
  
% Number of IPTG binding to LacI:
  n = p(17);
  

% ------------------------------------------------------------------------
% [2] MODEL - SYSTEM OF ODES
% ------------------------------------------------------------------------

% Ensuring concs do not become negative:
  if I <= 0
      I = 0;
  end
  if C <= 0
      C = 0;
  end
  
% IPTG in cell (I):
%   n = 2;
  sequestration = kf*(I^n)*L - kr*C;
  dx(1) = - n * sequestration;
  
% LacI expression and sequestration (dimer - L):
  dx(2) = - sequestration;
  
% LacI-IPTG-IPTG complex (C):
  dx(3) = sequestration;
