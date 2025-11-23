function g = forward_modeling(density_model, observation_points, block_size, G_constant)
% Forward Modeling for Rectangular Prism Gravity Anomaly
% Based on the analytical solution for rectangular prisms
% Reference: Camacho et al., Gravity inversion by means of growing bodies
% 
% Input:
%   density_model: 3D density array [nx, ny, nz] in kg/m^3
%   observation_points: Nx3 array of observation points [x, y, z]
%   block_size: [40, 40, 50] in meters
%   G_constant: Gravitational constant (6.674e-11 m^3 kg^-1 s^-2)
% 
% Output:
%   g: Calculated gravity anomaly at observation points [N]

% Parameters
if nargin < 4
    G_constant = 6.674e-11;
end

[nx, ny, nz] = size(density_model);
n_obs = size(observation_points, 1);
g = zeros(n_obs, 1);

% Block dimensions
dx = block_size(1);
dy = block_size(2);
dz = block_size(3);

% Loop over each block in the density model
for i = 1:nx
    for j = 1:ny
        for k = 1:nz
            % Get density of current block
            rho = density_model(i, j, k);
            
            if rho > 0  % Only process non-zero density blocks
                % Calculate block boundaries
                x_min = (i-1) * dx;
                x_max = i * dx;
                y_min = (j-1) * dy;
                y_max = j * dy;
                z_min = (k-1) * dz;
                z_max = k * dz;
                
                % Calculate gravity contribution from this block
                g = g + gravity_prism(observation_points, x_min, x_max, y_min, y_max, z_min, z_max, rho, G_constant);
            end
        end
    end
end

end

function g_prism = gravity_prism(obs_points, x1, x2, y1, y2, z1, z2, rho, G)
% Calculate gravity anomaly from a rectangular prism
% Using analytical formula

n_obs = size(obs_points, 1);
g_prism = zeros(n_obs, 1);

for p = 1:n_obs
    x0 = obs_points(p, 1);
    y0 = obs_points(p, 2);
    z0 = obs_points(p, 3);
    
    % Relative coordinates
    dx1 = x1 - x0;
    dx2 = x2 - x0;
    dy1 = y1 - y0;
    dy2 = y2 - y0;
    dz1 = z1 - z0;
    dz2 = z2 - z0;
    
    % Calculate potential field contribution
    g_val = 0;
    for i_sign = [-1, 1]
        for j_sign = [-1, 1]
            for k_sign = [-1, 1]
                dx = i_sign * (i_sign == -1 ? dx1 : dx2);
                dy = j_sign * (j_sign == -1 ? dy1 : dy2);
                dz = k_sign * (k_sign == -1 ? dz1 : dz2);
                
                r = sqrt(dx^2 + dy^2 + dz^2);
                
                if r > 1e-10
                    sign_factor = i_sign * j_sign * k_sign;
                    g_val = g_val + sign_factor * atan2(dx*dy, dz*r);
                end
            end
        end
    end
    
g_prism(p) = G * rho * g_val;
end

end