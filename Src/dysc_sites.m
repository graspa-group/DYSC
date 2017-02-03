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

classdef dysc_sites < handle
    
    properties (SetAccess = private) 
        coordinates=[];     %[double]   (Nx2) the matrix of spatial coordinates
        unit=[];            %[string]   (1x1) 'deg', 'km' or 'm'
        type='sparse';      %[string]   (1x1) 'sparse' or 'grid'
        grid_size=[];       %[integer]  (2x1) the number of rows and columns of the grid
        
        distance_matrix=[]; %[double]   (NxN) the distance matrix (in the unit of the spatial coordinates)
    end
    
    methods
        
        function obj = dysc_sites(coordinates,unit,type,grid_size)
            %DESCRIPTION: the constructor of the class dysc_sites
            %
            %INPUT
            %coordinates: [double] (Nx2) the matrix of spatial coordinates
            %unit:        [string] (1x1) 'deg', 'km' or 'm'
            %<type>:      [string] (1x1) 'sparse' or 'grid' (default: 'sparse')
            %<grid_size>: [integer](2x1) the grid size (default: [])
            %
            %OUTPUT
            %obj:         [dysc_sites object] (1x1) an object of class dysc_sites
            
            if nargin<2
                error('coordinates and unit must be provided');
            end
            
            if nargin==3
                if strcmpi(type,'grid')==1
                    error('grid_size must be provided');
                end
            end
            
            obj.coordinates=coordinates;
            obj.unit=unit;
            obj.type=type;
            
            if nargin>3
                obj.grid_size=grid_size;
            end
            
            obj.set_distance_matrix;
        end
        
        function set_distance_matrix(obj)
            if strcmpi(obj.unit,{'deg'})
                obj.distance_matrix=zeros(size(obj.coordinates,1));
                for i=1:size(obj.coordinates,1)
                    obj.distance_matrix(i,:)=distance(obj.coordinates(i,1),obj.coordinates(i,2),obj.coordinates(:,1),obj.coordinates(:,2));
                end
            end
            if strcmpi(obj.unit,{'km'})
                obj.distance_matrix=squareform(pdist(obj.coordinates,'euclidean'));
            end
            if strcmpi(obj.unit,{'m'})
                obj.distance_matrix=squareform(pdist(obj.coordinates,'euclidean'))/1000;
            end
        end

        function N = get_N(obj)
            N=size(obj.coordinates,1);
        end
        
        function set.coordinates(obj,coordinates)
            if not(size(coordinates,2)==2)
                error('coordinates must be a Nx2 matrix');
            end
            obj.coordinates=coordinates;
        end
        
        function set.unit(obj,unit)
            if sum(strcmpi(unit,{'deg','km','m'}))==0
                error('unit must be ''deg'' or ''km'' or ''m''');
            end
            obj.unit=unit;
        end
        
        function set.type(obj,type)
            if sum(strcmpi(type,{'sparse','grid'}))==0
                error('type must be ''sparse'' or ''grid''');
            end
            obj.type=type;
        end
        
        function set.grid_size(obj,grid_size)
            if not(length(grid_size)==2)
                error('grid_size must be a 2x1 vector');
            end
            if not(prod(grid_size)==size(obj.coordinates,1))
                error('The product of the elements of grid_size must be equal to the number of coordinates');
            end
            obj.grid_size=grid_size;
        end
    end
end
        