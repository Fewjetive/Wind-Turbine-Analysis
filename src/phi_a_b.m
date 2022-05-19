main()
function main()
    % [U_1, R, Omega, J, alpha, C_L, C_D] = read_data("NACA6409");
    % fprintf('U_1 = %f, R = %f, Omega = %f, J = %f, alpha = %f, C_L = %f, C_D = %f', U_1, R, Omega, J, alpha, C_L, C_D)
    [a, b, phi] = cal_a_b_phi(0.5, 0.5, "NACA6409");
    fprintf('(a = %f, b = %f, phi = %f)\n', a, b, phi);
end
function [U_1, R, Omega, J, alpha, C_L, C_D] = read_data(turbine)
    text = fileread("data\common_data.txt");
    text = split(text);
    i = 1;
    while ~strcmp(text(i, 1), "-")
        if text(i, 1) == "U_1"
            U_1 = str2double(text(i + 1, 1));
        elseif text(i, 1) == "R"
            R = str2double(text(i + 1, 1));
        elseif text(i, 1) == "Omega"
            Omega = str2double(text(i + 1, 1));
        elseif text(i, 1) == "J"
            J = str2double(text(i + 1, 1));
        end
        i = i + 2;
    end
    while ~strcmp(text(i, 1), turbine)
        i = i + 1;
    end
    i = i + 1;
    while 1
        if text(i, 1) == "alpha"
            alpha = str2double(text(i + 1, 1));
        elseif text(i, 1) == "C_L"
            C_L = str2double(text(i + 1, 1));
        elseif text(i, 1) == "C_D"
            C_D = str2double(text(i + 1, 1));
            break;
        end
        i = i + 2;
    end
end
function [a, b, phi] = cal_a_b_phi(r, sigma, turbine)
    if ~isfloat(r)
        error('r must be a double.')
    end
    if ~isfloat(sigma)
        error('sigma must be a double.')
    end
    if sigma < 0 || sigma > 1
        error('sigma must be in range [0, 1]. Input: %f', sigma)
    end
    [U_1, R, Omega, J, ~, C_L, C_D] = read_data(turbine); 
    % U_1 = 20; R = 1; Omega = 2 * pi; J = 1.256; C_L = 1.4553; C_D = 0.02363;
    D = 2 * R; a = 0.5; b = 0.5; i = 1;
    while true
        pa = a; pb = b;
        phi = cal_phi(a, b, r, U_1, Omega);
        a = cal_a(phi, sigma, pa, C_L, C_D);
        b = cal_b(phi, sigma, pa, r, C_L, C_D, J, D);
        i = i + 1;
        fprintf('(a = %f, b = %f, phi = %f)\n', a, b, phi);
        if i > 1000
            error('Claculate too many times!')
        elseif abs(a - pa) + abs(b - pb) <= 0.0000001
            break;
        end
    end
end
function phi = cal_phi(a, b, r, U_1, Omega)
    phi = atan(U_1 * (1 - a) / (Omega * r * (1 + b)));
end
function a = cal_a(phi, sigma, pa, C_L, C_D)
    a = sigma * (1 - pa) * (C_L * cos(phi) + C_D * sin(phi)) / (4 * (sin(phi)) ^ 2);
end
function b = cal_b(phi, sigma, pa, r, C_L, C_D, J, D)
    b = sigma * (1 - pa) * J * (C_L * sin(phi) - C_D * cos(phi)) / (4 * (sin(phi))^2 * 2 * pi * r / D);
end
