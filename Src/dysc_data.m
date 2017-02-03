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


classdef dysc_data < handle
  
    properties
        y=[];  %[double] (NxT) the observation matrix
        X=[];  %[double] (NxTxB) the covariate matrix
    end
    
    properties (SetAccess = private)
        N=[];  %[double] (1x1) the number of sites
        T=[];  %[double] (1x1) the number of temporal steps
        B=[];  %[double] (1x1) the number of covariates
    end
    
    methods
        
        function obj = dysc_data(y,X)
            %DESCRIPTION: is the constructor of the class dysc_data
            %
            %INPUT
            %y:         [double] (NxT) the observation matrix
            %X:         [double] (NxbxT) the covariate matrix
            %
            %OUTPUT
            %obj:       [dysc_data object] (1x1) an object of class dysc_data
            
            if nargin<1
                error('The input parameter y must be provided');
            end

            if nargin>=2
                if not(size(y,1)==size(X,1))
                    error('y and X must have the same number of rows');
                end
                if not(size(y,2)==size(X,2))
                    error('y and X must have the same number of temporal steps');
                end
            end
            N=size(y,1);
            T=size(y,2);

            obj.y=y;
            obj.N=N;
            obj.T=T;
            if nargin>=2
                obj.X=X;
                B=size(X,3);
                obj.B=B;
            else
                obj.B=0;
            end
        end
      end
end
        
        


