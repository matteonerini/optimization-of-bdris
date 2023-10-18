function [Theta, V] = func_theta(hRI_norm,hIT_norm,NG)
% Compute the scattering matrix Theta given the normalized channels
% hRI_norm and hIT_norm, and the group size NG
% Inputs:   hRI_norm: normalized channel hRI
%           hIT_norm: normalized channel hIT
%           NG: group size (NG=0 means fully connected)
% Outputs:  Theta: scattering matrix

NI = size(hIT_norm,1); % Number of RIS elements
if NG == 0 % Fully connected
    NG = NI;
end
G = NI/NG; % Number of groups

Theta = [];

if NG == 1 % Single connected

    % Compute Theta
    theta = - angle(hRI_norm) - angle(hIT_norm.');
    Theta = diag(exp(1i * theta));

else % Group or fully connected

    for g = 1:G
    
        % Truncated channels
        hRI_g = hRI_norm(NG*(g-1)+1:NG*g);
        hIT_g = hIT_norm(NG*(g-1)+1:NG*g);
        hRI_g_norm = hRI_g / norm(hRI_g);
        hIT_g_norm = hIT_g / norm(hIT_g);
    
        % Matrix A
        RRI = hRI_g_norm' * hRI_g_norm;
        RIT = hIT_g_norm * hIT_g_norm';
        ARI = (RRI + RRI.') / 2;
        AIT = (RIT + RIT.') / 2;
        A = ARI - AIT;
        
        % Eigenvalue decomposition of A
        [U,Delta] = eig(A);
        delta = flip(diag(Delta)); % Order delta in decreasing order
        
        % Compute matrix T distinguishing three cases
        T = zeros(NG);

        if NG == 2

            T = [sqrt(1/2), sqrt(1/2);
                 sqrt(1/2), -sqrt(1/2)];

        elseif NG == 3

            T = [sqrt(-delta(3)/(delta(1)-delta(3))), sqrt(delta(1)/(2*(delta(1)-delta(3)))), -sqrt(delta(1)/(2*(delta(1)-delta(3))));
                 0, sqrt(1/2), sqrt(1/2);
                 sqrt(delta(1)/(delta(1)-delta(3))), -sqrt(-delta(3)/(2*(delta(1)-delta(3)))), sqrt(-delta(3)/(2*(delta(1)-delta(3))))];

        else

            T(1,1) = sqrt(-delta(NG-1) / (delta(1) - delta(NG-1)));
            T(NG-1,1) = sqrt(delta(1) / (delta(1) - delta(NG-1)));
            T(2,2) = sqrt(-delta(NG) / (delta(2) - delta(NG)));
            T(NG,2) = sqrt(delta(2) / (delta(2) - delta(NG)));
            T(1,3) = sqrt(1/2) * T(NG-1,1);
            T(2,3) = sqrt(1/2) * T(NG,2);
            T(NG-1,3) = -sqrt(1/2) * T(1,1);
            T(NG,3) = -sqrt(1/2) * T(2,2);
            T(1,4) = sqrt(1/2) * T(NG-1,1);
            T(2,4) = -sqrt(1/2) * T(NG,2);
            T(NG-1,4) = -sqrt(1/2) * T(1,1);
            T(NG,4) = sqrt(1/2) * T(2,2);

            T(3:NG-2,5:NG) = eye(NG-4);

        end

        % Compute matrices V, D, and Theta
        V = flip(U,2) * T; % Order U according to delta
        hRI_g_bar = hRI_g_norm * V;
        hIT_g_bar = V' * hIT_g_norm;
        theta = - angle(hRI_g_bar) - angle(hIT_g_bar.');
        d = exp(1i * theta);
        D = diag(d);
        Theta_tmp = V * D * V';
        
        Theta = blkdiag(Theta,Theta_tmp); % The scattering matrix is block diagonal

    end

end

end