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

function [yout,tout,params,F] = multicell_LI(params)

% multicell_LI simulates lateral inhibition in a hexagonal lattice. 
% The structure params contains the model parameters of the system. 
% TOUT is a vector containing the time points of the solution 
% between 0 and Tmax. YOUT is a matrix containing the numerical 
% solution for each variable for each time point. Each row in 
% YOUT is a vector of the size of TOUT. F is a movie of the simulation. 

Tmax=30; tspan=[0 Tmax]; % set time for simulation

% get the default parameters if none provided
if(nargin < 1)
    params=defaultparams; 
end

P=params.P;     % number of cells per column
Q=params.Q;     % number of columns - MUST BE EVEN
k=P*Q;          % number of cells

% get the connectivity matrix
params.connectivity=getconnectivityM(P,Q);

% setting the initial conditions + noise
y0=getIC(params,k);

% run simulation with lateral inhibition
[tout,yout] = ode23(@li,tspan,y0,[],params);

% show time traces of two cells with lateral inhibition
plot2cells(tout,yout,k)
% 
% show lattice simulation
F=movielattice(tout,yout,P,Q,k);


function dy = li(t,y,params)

nu=params.nu;
betaD=params.betaD;
betaR=params.betaR;
h=params.h;
m=params.m;
M=params.connectivity;
k=length(M);

D = y(1:k); % levels of Delta in cells 1 to k
R = y(k+1:2*k); % levels of Repressor in cells 1 to k
Dneighbor=M*y(1:k);% average Delta level in the neighboring cells

% differential equations for Delta and repressor levels
dD = nu * (betaD.*1./(1 + R.^h)-D); 
dR = betaR.*Dneighbor.^m./(1 + Dneighbor.^m)-R; 
dy = [dD;dR]; 

function params=defaultparams

params.nu=1;        % ratio of degradation rates
params.betaD=50;    % normalized Delta production
params.betaR=50;    % normalized repressor production
params.h=3;         % Hill coefficient repression function
params.m=3;         % Hill coefficient activating function
params.sigma=0.2;% noise amplitude in initial conditions
params.P=18;        % number of cells per column
params.Q=18;        % number of columns - MUST BE EVEN

function M=getconnectivityM(P,Q)

k=P*Q;          % number of cells
M=zeros(k,k);   % connectivity matrix
w=1/6;          % weight for interactions

% calculating the connectivity matrix
for s=1:k
    kneighbor=findneighborhex(s,P,Q); 
    for r=1:6
        M(s,kneighbor(r))=w;
    end
end

function y0=getIC(params,k)

U=rand(k,1) - 1/2; % a uniform random distribution
epsilon=1e-5;   % multiplicative factor of Delta initial condition
D0=epsilon*params.betaD.*(1 + params.sigma*U); % initial Delta levels 
R0=zeros(k,1);  % initial repressor levels
y0=[D0;R0];     % vector of initial conditions

function plot2cells(tout,yout,k)

figure(21)
clf
for i=1:2
    subplot(1,2,i)
    plot(tout,yout(:,i),'-r','linewidth',2)   % plot Delta levels 
    hold on
    plot(tout,yout(:,k+i),'-b','linewidth',2) % plot repressor levels
    title(['cell #',num2str(i)])
    xlabel('t [a.u]'); ylabel('concentration [a.u]')
    legend('d','r')
end

function out = findneighborhex(ind,P,Q)

% This function finds the 6 neighbors of cell ind
[p,q] = ind2pq(ind,P);

% above and below:
out(1) = pq2ind(mod(p,P)+1,q,P);
out(2) = pq2ind(mod(p-2,P)+1,q,P);

% left and right sides:
qleft = mod(q-2,Q)+1;
qright = mod(q,Q)+1;

if q/2~=round(q/2),
    pup = p;
    pdown = mod(p-2,P)+1;
else 
    pup = mod(p,P)+1;
    pdown = p;
end;
out(3) = pq2ind(pup,qleft,P);
out(4) = pq2ind(pdown,qleft,P);
out(5) = pq2ind(pup,qright,P);
out(6) = pq2ind(pdown,qright,P);

function ind=pq2ind(p,q, P)
ind = p + (q-1)*P;

function [p,q]=ind2pq(ind, P)
q = 1+floor((ind-1)/P);
p = ind - (q-1)*P;

function plotHexagon(p0,q0,c)

% this function plots a hexagon centered at coordinates p,q

s32 = sqrt(3)/4;
q = q0*3/4;
p = p0*2*s32;
if q0/2 == round(q0/2),
   p = p+s32;
end;

x(1)=q-.5; x(2)=q-.25; x(3)=q+.25; 
x(4)=q+.5; x(5)=q+.25; x(6)=q-.25;

y(1)=p ; y(2)=p+s32; y(3)=p+s32; 
y(4)=p; y(5)=p-s32; y(6)=p-s32;

patch(x,y,c,'linewidth',2);

function F=movielattice(tout,yout,P,Q,k)

% This function generates a movie of patterning in hexagonal 
% lattice. The color represents the level of Delta. It also 
% saves the movie as an AVI file. 

figure(22)
Cmax=max(yout(end,1:k)); % finds max(Delta) at the end point
frameind=0;
for tind = 1:5:length(tout),   % shows every 5th frame
    clf;
    for i = 1:P,
        for j = 1:Q,
            ind = pq2ind(i,j,P);
            mycolor = min([yout(tind,ind)/Cmax,1]); 
            plotHexagon(i,j,[1-mycolor,1-mycolor,1]);
        end;
    end;
    axis image; axis off; box off;
    
    frameind=frameind+1;
    F(frameind) = getframe; % generates a movie variable
end;
% save movie in avi format
movie2avi(F,'movielattice','compression','none'); 
