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

function [yout,tout,params] = transcis2cell_LI(params)

% transcis2cell_LI simulates trans-annihilation with cis-inactivation
% between two cells. The structure params contains the model parameters of the system. 
% TOUT is a vector containing the time points of the solution 
% between 0 and Tmax. YOUT is a matrix containing the numerical 
% solution for each variable for each time point. Each row in 
% YOUT is a vector of the size of TOUT.

Tmax=100; tspan=[0 Tmax]; % set time for simulation
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
betaN=params.betaN;
betaR=params.betaR;
m=params.m;
h=params.h;
M=params.connectivity;
k=length(M);
mu=params.mu;
kc=params.kc;
kt=params.kt;

D = y(1:k);         % levels of Delta in cells 1 to k
R = y(k+1:2*k);     % levels of Repressor in cells 1 to k
N = y(2*k+1:3*k);   % levels of Repressor in cells 1 to k
Dneighbor=M*y(1:k);       % Delta level in the neighboring cells
Nneighbor=M*y(2*k+1:3*k); % Notch level in the neighboring cells

% differential equations for Delta, repressor, and Notch levels
dN = mu * (betaN - kt.*N.*Dneighbor-kc.*N.*D-N); 
dD = nu * (betaD.*1./(1 + R.^h)-kt.*D.*Nneighbor-kc.*N.*D-D); 
dR = betaR.*(kt.*N.*Dneighbor).^m./(1 + (kt.*N.*Dneighbor).^m)-R; 
dy = [dD;dR;dN];  

function params=defaultparams

params.nu=1;
params.betaD=50;
params.betaN=1;
params.betaR=200;
params.m=1;
params.h=1;
params.sigma=0.2;
params.mu=1;
params.kc=10;
params.kt=1;

function M=getconnectivityM

M=[0 1;1 0]; % 2 cell connectivity matrix

function y0=getIC(params,k)

U=rand(k,1) - 1/2; % a uniform random distribution
epsilon=1e-5;   % multiplicative factor of Delta initial condition
D0=epsilon*params.betaD.*(1 + params.sigma*U); % initial Delta levels 
R0=zeros(k,1);  % initial repressor levels
N0=params.betaN.*ones(k,1);  % initial Notch levels are betaN
y0=[D0;R0;N0];  % vector of initial conditions

function plot2cells(tout,yout,k)

figure(21)
clf
for i=1:2
    subplot(1,2,i)
    plot(tout,yout(:,i),'-r','linewidth',2)   % plot Delta levels 
    hold on
    plot(tout,yout(:,k+i),'-b','linewidth',2) % plot repressor levels
    title(['cell #',num2str(i)])
    xlabel('time [a.u]'); ylabel('concentration [a.u]')
    legend('d','r')
end
