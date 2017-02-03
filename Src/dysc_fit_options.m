%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DYSC - DYnamic Spatiotemporal Clustering                                  %
%%%                                                                           %
%%% Authors: Francesco Finazzi and Lucia Paci                                 %
%%% E-mail: francesco.finazzi@unibg.it - lucia.paci@unicatt.it                %
%%% Affiliation: University of Bergamo - Università  Cattolica del Sacro Cuore%
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


classdef dysc_fit_options < handle
    
    properties
        iterations=10000;         %[integer>0]    (1x1) the number of MCMC iterations
        burn_in=2500;             %[integer>=0]   (1x1) the number of burn-in iterations
        thin=1;                   %[integer>=0]   (1x1) thinning
        
        flag_update_sigma=1;      %[boolean]      (1x1) 0: the parameter is not estimated; 1:the parameter is estimated;
        flag_update_lambda=1;     %[boolean]      (1x1) 0: the parameter is not estimated; 1:the parameter is estimated;
        flag_update_tau=1;        %[boolean]      (1x1) 0: the parameter is not estimated; 1:the parameter is estimated;
        flag_update_z=1;          %[boolean]      (1x1) 0: the parameter is not estimated; 1:the parameter is estimated;
        flag_update_G=1;          %[boolean]      (1x1) 0: the parameter is not estimated; 1:the parameter is estimated;
        flag_update_rho=1;        %[boolean]      (1x1) 0: the parameter is not estimated; 1:the parameter is estimated;
        flag_update_phi=1;        %[boolean]      (1x1) 0: the parameter is not estimated; 1:the parameter is estimated;
        flag_update_w=1;          %[boolean]      (1x1) 0: the parameter is not estimated; 1:the parameter is estimated;  
        flag_update_beta=1;       %[boolean]      (1x1) 0: the parameter is not estimated; 1:the parameter is estimated;
        flag_update_theta_zeta=0; %[boolean]      (1x1) 0: the parameter is not estimated; 1:the parameter is estimated;
        
        flag_relabelling=1;       %[boolean]      (1x1) 0: the relabelling of z is off; 1: the relabelling of z is on
        
        proposal_type='rw';       %[string]       (1x1) 'rw': random-walk proposal, 'conditiona': conditional proposal
        
        tuning_phi=0.01;          %[double]       (K-1x1) the tuning parameter for the proposal on phi. If a scalar is provided, the vector K-1x1 is created using the scalar
        tuning_beta=0.01;         %[double]       (Bx1) the tuning parameter for the proposal on beta. If a scalar is provided, the vector Bx1 is created using the scalar
    end
    
    methods
        
        function obj = dysc_fit_options(iterations,burn_in,thin)
            %DESCRIPTION: is the constructor of the class dysc_fit_options
            %
            %INPUT
            %<iterations>:  [integer>0]    (1x1) the number of MCMC iterations
            %<burn_in>:     [integer>=0]   (1x1) the number of burn-in iterations
            %<thin>:        [integer>=0]   (1x1) thinning
            %
            %OUTPUT
            %obj:           [dysc_fit_options object] (1x1) an object of class dysc_fit_options
            
            if nargin>=1
                obj.iterations=iterations;
            end
            %burn in
            if nargin>=2
                if burn_in>obj.iterations
                    error('burn_in must be lower than iterations');
                end
                obj.burn_in=burn_in;
            end
            %thinning
            obj.thin=thin;
        end
        
                
        function set.iterations(obj,iterations)
            if iterations<1
                error('iterations must be > 0');
            end
            obj.iterations=iterations;
        end
        
        function set.burn_in(obj,burn_in)
            if burn_in<1
                error('burn_in must be > 0');
            end
            if obj.iterations<burn_in
                error('burn_in must be lower than iterations');
            end
            obj.burn_in=burn_in;
        end
        
            function set.thin(obj,thin)
            if thin<1
                error('thinning must be > 0');
            end
            obj.thin=thin;
        end
        
        function set.flag_update_sigma(obj,flag_update_sigma)
            if not(flag_update_sigma==0||flag_update_sigma==1)
                error('flag_update_sigma must be either 0 or 1');
            end
            obj.flag_update_sigma=flag_update_sigma;
        end
        
        function set.flag_update_lambda(obj,flag_update_lambda)
            if not(flag_update_lambda==0||flag_update_lambda==1)
                error('flag_update_lambda must be either 0 or 1');
            end
            obj.flag_update_lambda=flag_update_lambda;
        end
        
        function set.flag_update_tau(obj,flag_update_tau)
            if not(flag_update_tau==0||flag_update_tau==1)
                error('flag_update_tau must be either 0 or 1');
            end
            obj.flag_update_tau=flag_update_tau;
        end
        
        function set.flag_update_z(obj,flag_update_z)
            if not(flag_update_z==0||flag_update_z==1)
                error('flag_update_z must be either 0 or 1');
            end
            obj.flag_update_z=flag_update_z;
        end
        
        function set.flag_update_G(obj,flag_update_G)
            if not(flag_update_G==0||flag_update_G==1)
                error('flag_update_sigma must be either 0 or 1');
            end
            obj.flag_update_G=flag_update_G;
        end
        
        function set.flag_update_rho(obj,flag_update_rho)
            if not(flag_update_rho==0||flag_update_rho==1)
                error('flag_update_rho must be either 0 or 1');
            end
            obj.flag_update_rho=flag_update_rho;
        end
        
        function set.flag_update_phi(obj,flag_update_phi)
            if not(flag_update_phi==0||flag_update_phi==1)
                error('flag_update_phi must be either 0 or 1');
            end
            obj.flag_update_phi=flag_update_phi;
        end
        
        function set.flag_update_w(obj,flag_update_w)
            if not(flag_update_w==0||flag_update_w==1)
                error('flag_update_w must be either 0 or 1');
            end
            obj.flag_update_w=flag_update_w;
        end
        
        function set.flag_update_beta(obj,flag_update_beta)
            if not(flag_update_beta==0||flag_update_beta==1)
                error('flag_update_sigma must be either 0 or 1');
            end
            obj.flag_update_beta=flag_update_beta;
        end
        
        function set.flag_update_theta_zeta(obj,flag_update_theta_zeta)
            if not(flag_update_theta_zeta==0||flag_update_theta_zeta==1)
                error('flag_update_theta_zeta must be either 0 or 1');
            end
            obj.flag_update_theta_zeta=flag_update_theta_zeta;
        end  
       
        function set.flag_relabelling(obj,flag_relabelling)
            if not(flag_relabelling==0||flag_relabelling==1)
                error('flag_relabelling must be either 0 or 1');
            end
            obj.flag_relabelling=flag_relabelling;
        end
        
        function set.proposal_type(obj,proposal_type)
            if sum(strcmp(proposal_type,{'rw','conditional'}))==0
                error('proposal_type must be ''rw'' or ''conditional''');
            end
            obj.proposal_type=proposal_type;
        end
        
        function set.tuning_phi(obj,tuning_phi)
            if not(tuning_phi>0)
                error('tuning_phi must be > 0');
            end
            obj.tuning_phi=tuning_phi;
        end
        
        function set.tuning_beta(obj,tuning_beta)
            if not(tuning_beta>0)
                error('tuning_beta must be > 0');
            end
            obj.tuning_beta=tuning_beta;
        end
        
    end
end