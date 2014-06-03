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

function [yout,tout,params,F] = transcis_multicell_LI(params)

% transcis_multicell_LI simulates trans-annihilation with cis-inactivation 
% in a hexagonal lattice. The structure params contains the model parameters 
% of the system. TOUT is a vector containing the time points of the solution 
% between 0 and Tmax. YOUT is a matrix containing the numerical solution for 
% each variable for each time point. Each row in YOUT is a vector of the 
% size of TOUT. F is a movie showing the simulation. 

Tmax=100; tspan=[0 Tmax]; % set time for simulation

if(nargin < 1)
    params=defaultparams; % get the default parameters if none provided
end

P=params.P;  % number of cells per column
Q=params.Q;  % number of columns - MUST BE EVEN
k=P*Q; % number of cells

% get the connectivity matrix
params.connectivity=getconnectivityM(P,Q);

% setting the initial conditions + noise
y0=getIC(params,k);

% run simulation with lateral inhibition
[tout,yout] = ode15s(@li,tspan,y0,[],params);

% show time traces of two cells with lateral inhibition

plot2cells(tout,yout,k)

% show lattice simulation

F=movielattice(tout,yout,P,Q,k);

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

D = y(1:k); % levels of Delta in cells 1 to k
R = y(k+1:2*k); % levels of Repressor in cells 1 to k
N = y(2*k+1:3*k); % levels of Repressor in cells 1 to k
Dneighbor=M*y(1:k);% average Delta level in the neighboring cells
Nneighbor=M*y(2*k+1:3*k); % average Notch level in the neighboring cells

dN = mu * (betaN - kt.*N.*Dneighbor-kc.*N.*D-N); % differential equation describing Notch levels
dD = nu * (betaD.*1./(1 + R.^h)-kt.*D.*Nneighbor-kc.*N.*D-D); % differential equation describing Delta levels
dR = betaR.*(N.*Dneighbor).^m./(1 + (N.*Dneighbor).^m)-R; % differential equation describing repressor levels

dy = [dD;dR;dN];  

function params=defaultparams

params.nu=1;
params.betaD=50;
params.betaN=1;
params.betaR=200;
params.m=1;
params.h=1;
params.sigma=0.2;
params.P=12;
params.Q=12;
params.mu=1;
params.kc=10;
params.kt=1;

function M=getconnectivityM(P,Q)

k=P*Q; % number of cells
M=zeros(k,k); % connectivity matrix
w=1/6; % weight for interactions

% calculating the connectivity matrix
for s=1:k
    kneighbor=findneighborhex(s,P,Q); % finds the neighbors of cell s
    for r=1:6
        M(s,kneighbor(r))=w;
    end
end

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
    plot(tout,yout(:,i),'-r','linewidth',2)   % plot D levels 
    hold on
    plot(tout,yout(:,k+i),'-b','linewidth',2) % plot R levels
    title(['cell #',num2str(i)])
    xlabel('time [a.u]'); ylabel('concentration [a.u]')
    legend('d','r')
end

function out = findneighborhex(ind,P,Q)
[p,q] = ind2pq(ind,P);

%above and below:
out(1) = pq2ind(mod(p,P)+1,q,P);
out(2) = pq2ind(mod(p-2,P)+1,q,P);

%left side:
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

% This function plots a hexagon centered at hex lattice coordinates p,q

s32 = sqrt(3)/4;

q = q0*3/4;
p = p0*2*s32;

if q0/2 == round(q0/2),
   p = p+s32;
end;

x(1) = q-.5; x(2) = q-.25; x(3) = q+.25; x(4) = q+.5; x(5) = q+.25; x(6) = q-.25;

y(1) = p ; y(2) = p+s32; y(3) = p+s32; y(4) = p; y(5) = p-s32; y(6) = p-s32;

c=min(c,ones(1,3));

patch(x,y,c,'linewidth',2);

function F=movielattice(tout,yout,P,Q,k)

figure(22)
sy1 = sort(yout(end,1:k));
Norm = sy1(round(length(sy1)*0.95)); % find the Delta level in the high Delta cells
frameind=0;
for tind = 1:5:length(tout),   % shows every 5th frame
    clf;
    for i = 1:P,
        for j = 1:Q,
            ind = pq2ind(i,j,P);
            mycolor = min([yout(tind,ind)/ Norm,1]); % defined the normalized color of cell
            plotHexagon(i,j,[1-mycolor,1-mycolor,1]);
        end;
    end;
    axis image; axis off; box off;
    
    frameind=frameind+1;
    F(frameind) = getframe; % generates a movie variable
end;
movie2avi(F,'movielattice'); % save movie in avi format
