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

function paramsearch_LI

% This code plots the log(Dmax/Dmin) and the time required for patterning 
% for different betaD and betaR. TO RUN FASTER, COMMENT OUT PLOT2CELLS AND 
% MOVIELATTICE within multicell_LI function. For a better exploration in 
% the parameter space, it is important to set a large enough Tmax parameter 
% in multicell_LI function.


% fixed parameters
params.nu=1;
params.h=3;
params.m=3;
params.sigma=0.2;
params.Q=12;
params.P=12;
k=params.Q*params.P;

% variable parameters
betaD=logspace(-1,3,20); % creates a series of betaD from 0.1 to 1000
betaR=logspace(-1,3,20); % creates a series of betaR from 0.1 to 1000

ind=0; h=waitbar(0,'% of progress'); % generates a waitbar
for i=1:length(betaD)
    for j=1:length(betaR)
        params.betaD=betaD(i);
        params.betaR=betaR(j);
        ind=ind+1; waitbar(ind/(length(betaD)*length(betaR))) 
        [yout,tout] = multicell_LI(params); % calling the LI solver
        
        % finding max and min values of D
        Dmax(i,j)=max(max(yout(end,1:k)));
        Dmin(i,j)=abs(min(min(yout(end,1:k))));
        
        % finding cases where patterning occurs (when Dmax/Dmin>1.2)
        % and getting the patterning time 
        if Dmax(i,j)/Dmin(i,j)>1.2 
            T(i,j)=getPatterningTime(tout,yout,k,Dmax(i,j),Dmin(i,j));
        else
            T(i,j)=NaN; % patterning time is not set for the no patterning case
        end
    end
end
close(h)
figure(23)
imagesc(log10(betaD),log10(betaR),log10(Dmax./Dmin));
set(gca,'YDir','normal')
xlabel('log(\beta_d)','fontsize',14);
ylabel('log(\beta_r)','fontsize',14);
title('log(d_{max}/d_{min})','fontsize',14)
colorbar

figure(24)
imagesc(log10(betaD),log10(betaR),T);
set(gca,'YDir','normal')
xlabel('log(\beta_d)','fontsize',14);
ylabel('log(\beta_r)','fontsize',14);
title('Time for patterning [normalized time]','fontsize',14)
colorbar

function T=getPatterningTime(tout,yout,k,Dmax,Dmin)

% This function estimates patterning time. This is done by the follwoing 3 
% steps:
% 1. find all the high D cells ('onCells')
% 2. find the time it takes for each 'on cell' to reach 90% of its final
% level ('TonCells')
% 3. get median value of the times calculated in stage 2

onCells=find(yout(end,1:k)>0.5*(Dmax+Dmin));
for i=1:length(onCells)
    Tind=find(yout(:,onCells(i))>0.9*yout(end,onCells(i)),1,'first');
    TonCells(i)=tout(Tind);
end 

T=median(TonCells);
