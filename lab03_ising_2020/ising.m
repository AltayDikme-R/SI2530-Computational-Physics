
% Monte Carlo simulation of the 2D ISing model

%parameters
L=30                 % size = L*L
sweeps=400             % number of update sweeps through the lattice
Tc=2/log(1+sqrt(2)); % critical temperature
T=0.5*Tc             % temperature
%T=3*Tc               % temperature
h=0                  % applied magnetic field
%h=1                  % applied magnetic field


% initialize spin S, magnetization M, and energy E
clear S;
for i=1:L
for j=1:L
   S(i,j) = 1;
%  random initial configuration
   if rand>0.5
      S(i,j)=-1;
   end
end
end
% These formulas only work for all spins up:
%M=L*L;
%E=-2*L*L-h*M;


% index tables for nearest neighbors using periodic boundary
% conditions
for i=1:L
   neigm(i)=i-1;
   neigp(i)=i+1;
end
neigm(1)=L;
neigp(L)=1;


% initialize result vectors and time step
%clear temp;
%clear energy;
%clear heatcap;
%clear mag;
%clear absmag;
%clear mag2;
%clear mag4;
%step=0;


% loop over temperatures
%for T = 1.5 : 0.1 : 4
%T


% Initialize averages to zero
%energy_ave=0;
%energy2_ave=0;
%mag_ave=0;
%absmag_ave=0;
%mag2_ave=0;
%mag4_ave=0;


% MC sweep through the lattice
for isweeps=1:sweeps
%   m = 0;
   slump = rand(L,L);
   for it=1:L*L
      i=fix(L*rand)+1;
      j=fix(L*rand)+1;
%   for i=1:L
%   for j=1:L
      dE = 2*S(i,j)*(S(neigm(i),j) + S(neigp(i),j) + S(i,neigm(j)) + S(i,neigp(j)) + h);
      w = exp(-dE/T);
      if w >= slump(i,j)
         S(i,j)=-S(i,j);
%         M=M+2*S(i,j);
%         E=E+dE;
      end
%   end
    end

% After initial sweeps have been discarded, accumulate to averages
%   if isweeps>sweeps/10
%      energy_ave=energy_ave+E;
%      energy2_ave=energy2_ave+E^2;
%      mag_ave=mag_ave+M;
%      absmag_ave=absmag_ave+abs(M);
%      mag2_ave=mag2_ave+M^2;
%      mag4_ave=mag4_ave+M^4;
%   end

% plot state
   pcolor(S)
   axis equal
   pause(0.01)
end


% compute averages and store in vectors for plotting
%step=step+1;
%temp(step) = T;
%energy_ave=energy_ave/(sweeps*9/10);
%energy2_ave=energy2_ave/(sweeps*9/10);
%heatcap(step)=(energy2_ave-energy_ave^2)/T^2/L/L;
%energy(step)=energy_ave/L/L;
%mag(step)=mag_ave/(sweeps*9/10)/L/L;
%absmag(step)=absmag_ave/(sweeps*9/10)/L/L;
%mag2_ave=mag2_ave/(sweeps*9/10);
%mag4_ave=mag4_ave/(sweeps*9/10);
%binder(step)=1-mag4_ave/mag2_ave^2/3;
%end


% plot magnetization
%plot(temp,mag,'r-',temp,mag,'bo')
%title('size L=10, magnetic field h=0.5')
%xlabel('temperature T')
%ylabel('magnetization m')








