%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DYSC - DYnamic Spatiotemporal Clustering                                  %
%%%                                                                           %
%%% Authors: Francesco Finazzi and Lucia Paci                                 %
%%% E-mail: francesco.finazzi@unibg.it - lucia.paci@unicatt.it                %
%%% Affiliation: University of Bergamo - Università Cattolica del Sacro Cuore %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This file is part of DYSC.
%
% DYSC is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 2 of the License, or
% (at your option) any later version.
%
% DYSC is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with DYSC. If not, see <http://www.gnu.org/licenses/>.


classdef dysc_par < handle
    
    properties
        K         = 2;   %[integer>0]   (1x1) the number of clusters (default=2)
        
        sigma     =[];   %[double]      (1x1)
        G         =[];   %[double]      (Kx1)
        tau       =[];   %[double]      (Kx1)
        rho       =[];   %[double]      (K-1x1) 
        lambda    =[];   %[double]      (K-1x1) 
        z         =[];   %[double]      (KxT+1)
        phi       =[];   %[double]      (N*KxT+1) 
        w         =[];   %[double]      (NxT+1) 
        beta      =[];   %[double]      (BxK)
        theta_zeta=[];   %[double>0]    (1x1)
    end
    
    properties (SetAccess = private, GetAccess = private)
        N=[];       %[integer>0]  (1x1) 
        T=[];       %[integer>0]  (1x1)
        B=[];       %[integer>=0] (1x1)
    end
    
    methods
        
        function obj = dysc_par(dysc_data,K)
            %DESCRIPTION: is the constructor of the class dysc_par
            %
            %INPUT
            %dysc_data:     [dysc_data object]  (1x1) an object of class dysc_data
            %
            %OUTPUT
            %obj:           [dysc_par object]   (1x1) an object of class dysc_par
            if nargin<1
                error('The input argument dysc_data must be provided');
            end
            
            if nargin>1
                obj.K=K;
            end
            
            N=dysc_data.N;
            obj.N=N;
            T=dysc_data.T;
            obj.T=T;
            B=dysc_data.B;
            obj.B=B;
            
            obj.sigma=1;
            obj.G=ones(obj.K,1);
            obj.tau=ones(obj.K,1);
            obj.z=zeros(obj.K,T+1);
            obj.rho=ones(obj.K-1,1);
            obj.lambda=ones(obj.K-1,1);
            obj.phi=zeros(N*obj.K,T+1);
            obj.w=ones(N,T+1);
            if B>0
                obj.beta=zeros(B,obj.K);
            end
        end
        
        function set.K(obj,K)
            if not(K==round(K))
                error('K must be integer');
            end
            if K<2
                error('K must be >= 2');
            end
            obj.K=K;
        end
        
        function set.sigma(obj,sigma)
            if not(size(sigma,1)==1)
                error('sigma must be a scalar');
            end
            if sum(sigma<=0)>0
                error('sigma must be > 0');
            end
            obj.sigma=sigma;
        end
        
        function set.G(obj,G)
            if not(size(G,1)==obj.K&&size(G,2)==1)
                error('G must be Kx1');
            end
            obj.G=G;
        end
        
        function set.tau(obj,tau)
            if not(iscell(tau))
                if sum(tau<=0)>0
                    error('All the elements of tau must be > 0');
                end
            else
                if not(length(tau))==obj.K
                    error('The number of cells in tau must be K');
                end
                for i=1:length(tau)
                    if sum(tau{i}<=0)>0
                        error('The elements in each cell of tau must be all > 0');
                    end
                end
            end
            obj.tau=tau;
        end
        
        function set.rho(obj,rho)
            if not(size(rho,1)==obj.K-1&&size(rho,2)==1)
                error('rho must be K-1x1');
            end
            obj.rho=rho;
        end
        
        function set.lambda(obj,lambda)
            if not(size(lambda,1)==obj.K-1&&size(lambda,2)==1)
                error('lambda must be K-1x1');
            end
            if sum(lambda<=0)>0
                error('All the elements of lambda must be > 0');
            end
            obj.lambda=lambda;
        end
        
        function set.z(obj,z)
            if not(size(z,1)==obj.K&&size(z,2)==(obj.T+1))
                error('z must be KxT+1');
            end
            obj.z=z;
        end
        
        function set.phi(obj,phi)
            if not(size(phi,1)==obj.N*obj.K&&size(phi,2)==obj.T+1)
                error('phi must be N*KxT+1');
            end
            obj.phi=phi;
        end    
        
        function set.w(obj,w)
            if not(size(w,1)==obj.N&&size(w,2)==obj.T+1)
                error('w must be NxT+1');
            end
            obj.w=w;
        end
        
        function set.beta(obj,beta)
            if not(size(beta,1)==obj.B&&size(beta,2)==obj.K)
                error('beta must be BxK');
            end
            obj.beta=beta;
        end
        
        function set.theta_zeta(obj,theta_zeta)
            if theta_zeta<=0
                error('theta_zeta must be > 0');
            end
            obj.theta_zeta=theta_zeta;
        end
        
    
    end
end