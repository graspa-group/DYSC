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


classdef dysc_hyper < handle
  
    properties
        %nu0 e S0 per Delta
        %prior di phi0 a priori normale multivariata media mu_phi0=0 e sigma_phi0=1 che poi moltiplicherà sigma_phi
        
        a_sigma=2;          %[double]      (1x1) the a parameter of the inverse gamma prior distribution on sigma
        b_sigma=1;          %[double]      (1x1) the b parameter of the inverse gamma prior distribution on sigma
        a_lambda=2;         %[double]      (1x1) the a parameter of the inverse gamma prior distribution on lambda
        b_lambda=1;         %[double]      (1x1) the b parameter of the inverse gamma prior distribution on lambda
        a_tau=2;            %[double]      (1x1) the a parameter of the inverse gamma prior distribution on tau
        b_tau=1;            %[double]      (1x1) the b parameter of the inverse gamma prior distribution on tau
        mu_z0=0;            %[double]      (Kx1) the mean of the normal prior distribution on z(0). If a scalar is provided, a vector Kx1 is built using the scalar
        sigma_z0=1000;      %[double>0]    (KxK) the variance of the normal prior distribution on z(0). If a scalar is provided, a diagonal matrix KxK is built using the scalar
        mu_beta0=0;         %[double]      (1x1) the mean of the normal prior distributions on beta
        sigma_beta0=1000;   %[double>0]    (1x1) the variance of the normal prior distributions on beta
        mu_phi0=0;          %[double]      (1x1) the mean of the normal prior distributions on beta
        sigma_phi0=1;       %[double>0]    (1x1) the variance of the normal prior distributions on beta
        
        sigma_G=10^4;       %[double>0]    (1x1) the variance of the prior on G
        sigma_rho=10^4;     %[double>0]    (1x1) the variance of the prior on rho
    end
    
    methods
        
        function obj = dysc_hyper()
            %DESCRIPTION: is the constructor of the class dysc_hyper
            %
            %OUTPUT
            %obj:      [dysc_hyper object] (1x1) an object of class dysc_hyper
        end

        function set.a_sigma(obj,a_sigma)
            if a_sigma<=0
                error('a_sigma must be > 0');
            end
            obj.a_sigma=a_sigma;
        end
        
        function set.a_lambda(obj,a_lambda)
            if a_lambda<=0
                error('a_lambda must be > 0');
            end
            obj.a_lambda=a_lambda;
        end
        
        function set.a_tau(obj,a_tau)
            if a_tau<=0
                error('a_tau must be > 0');
            end
            obj.a_tau=a_tau;
        end
        
        function set.b_sigma(obj,b_sigma)
            if b_sigma<=0
                error('b_sigma must be > 0');
            end
            obj.b_sigma=b_sigma;
        end
        
        function set.b_lambda(obj,b_lambda)
            if b_lambda<=0
                error('b_lambda must be > 0');
            end
            obj.b_lambda=b_lambda;
        end
        
        function set.b_tau(obj,b_tau)
            if b_tau<=0
                error('b_tau must be > 0');
            end
            obj.b_tau=b_tau;
        end
        
        function set.sigma_z0(obj,sigma_z0)
            if sum(diag(sigma_z0)<=0)>0
                error('sigma_z0 must be definite positive');
            end
            obj.sigma_z0=sigma_z0;
        end
        
        function set.sigma_beta0(obj,sigma_beta0)
            if sigma_beta0<=0
                error('sigma_beta0 must be > 0');
            end
            obj.sigma_beta0=sigma_beta0;
        end
        
        function set.sigma_phi0(obj,sigma_phi0)
            if sigma_phi0<=0
                error('sigma_phi0 must be > 0');
            end
            obj.sigma_phi0=sigma_phi0;
        end

        function set.sigma_rho(obj,sigma_rho)
            if sigma_rho<=0
                error('sigma_rho must be > 0');
            end
            obj.sigma_rho=sigma_rho;
        end

        function set.sigma_G(obj,sigma_G)
            if sigma_G<=0
                error('sigma_G must be > 0');
            end
            obj.sigma_G=sigma_G;
        end
        
      end
end
       