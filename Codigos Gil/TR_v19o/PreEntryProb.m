%//////////////////////////////////////////////////////////////////////////
% Source of this file: Ulyssea (2018), AER
%//////////////////////////////////////////////////////////////////////////


function [Prob]=PreEntryProb(q,sigma,dS,State)

% This function uses a method very close to Tauchen's (1986) to compute the
% pre-entry expected value in both sectors. Each point of the productivity
% space has a correspondant vector of probabilities.
%
% The method uses equispaced grid points instead of equal areas over the the grid (as presented in
% Adda's and Cooper's book, Appendix to Ch.3).
%
% Syntax:
%
% q           = pre-entry productivity of individual i.
% sigma       = standard deviation of the unantecipated shock.
% dS          = distance between grid points.
% State       = vector with the grid values of the discretized state space. 

S=length(State);
Prob=zeros(1,S);


% The end points:
x = (State(1,1)-q + dS/2)/sigma;
Prob(1,1)=normcdf(x,0,1);

x = (State(S,1)-q - dS/2)/sigma;
Prob(1,S) = 1 - normcdf(x,0,1);

% The remaining of the vector:

for j=2:1:S-1
    
    x1=(State(j,1)-q - dS/2)/sigma;
    x2=(State(j,1)-q + dS/2)/sigma;
        
    Prob(1,j) = normcdf(x2,0,1) - normcdf(x1,0,1);
    
end