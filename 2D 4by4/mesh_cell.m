classdef mesh_cell
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        
        rays struct
        scaler_flux
    end

    methods
        function obj = mesh_cell()
            
          
           obj.rays=struct("start_position",[],"end_position",[],"length",[],"angle",[],"angular_flux",[]);
           

        end

    end
end