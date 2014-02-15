% Copyright (C) 2014, David Sprinzak
% This program is part of Lateral Inhibition Tutorial.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, version 3 of the License.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [yout,tout,params] = twocell_LI(params)

% twocell_LI simulates lateral inhibition between two cells. The 
% structure params contains the model parameters of the system. 
% TOUT is a vector containing the time points of the solution 
% between 0 and Tmax. YOUT is a matrix containing the numerical 
% solution for each variable for each time point. Each row in 
% YOUT is a vector of the size of TOUT.  

Tmax=40; tspan=[0 Tmax]; % set time for simulation
k=2; % number of cells

% get the default parameters if none provided
if(nargin < 1)
    params=defaultparams;  
end

% get the connectivity matrix
params.connectivity=getconnectivityM;

% setting the initial conditions + noise
y0=getIC(params,k);

% run simulation with lateral inhibition
[tout,yout] = ode23(@li,tspan,y0,[],params);

% show time traces of two cells with lateral inhibition
plot2cells(tout,yout,k)


function dy = li(t,y,params) 

nu=params.nu;
betaD=params.betaD;
betaR=params.betaR;
h=params.h;
m=params.m;
M=params.connectivity;
k=length(M);

D = y(1:k);     % levels of Delta in cells 1 to k
R = y(k+1:2*k); % levels of repressor in cells 1 to k
Dneighbor=M*y(1:k);% Delta level in the neighboring cells

%differential equations for Delta and repressor levels
dD = nu * (betaD.*1./(1 + R.^h)-D); 
dR = betaR.*Dneighbor.^m./(1 + Dneighbor.^m)-R;
dy = [dD;dR];  

function params=defaultparams
 
params.nu=1;        % ratio of degradation rates
params.betaD=50;    % non-dimensional  Delta production
params.betaR=50;    % non-dimensional repressor production
params.h=3;         % Hill coefficient repression function
params.m=3;         % Hill coefficient activating function
params.sigma=0.2;   % noise amplitude in initial conditions
 
function M=getconnectivityM
 
M=[0 1;1 0]; % 2 cell connectivity matrix
 
function y0=getIC(params,k)
 
U=rand(k,1) - 1/2; % a uniform random distribution
epsilon=1e-5;   % multiplicative factor of Delta initial condition
D0=epsilon*params.betaD.*(1 + params.sigma*U); % initial Delta levels 
R0=zeros(k,1);  % initial repressor levels
y0=[D0;R0];     % vector of initial conditions

function plot2cells(tout,yout,k)

figure(21); clf
for i=1:2
    subplot(1,2,i)
    plot(tout,yout(:,i),'-r','linewidth',2)   %plot D levels 
    hold on
    plot(tout,yout(:,k+i),'-b','linewidth',2) %plot R levels
    title(['cell #',num2str(i)])
    xlabel('time [a.u]'); ylabel('concentration [a.u]')
    legend('d','r')
end
