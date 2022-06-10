main()
file_name = 'parameters.xlsx';
function main()
    r_num = 100;
    sigma_num = 1000;
    r_t = 14; % [1, r_t-1], SG6043, [r_t, r_num] -> NACA6409
    [B, U_1, R, D, Omega, rho, n, J, SIGMA] = read_data();
    [theta, naca_matrix, sg_matrix, NACA6409, SG6043, ~, ~] = get_standard(r_num, sigma_num, r_t, R, U_1, Omega, J, rho, B, D, n, SIGMA);   
    [C_T, C_Q, T, Q] = vary(r_num, sigma_num, r_t, R, SIGMA, U_1, Omega, 0.74, rho, B, D, n, theta, naca_matrix, sg_matrix, NACA6409, SG6043)
    % start = 0.452;
    % finish = 2.262;
    % % num = 100;
    % % increment = (finish - start) / num;
    % increment = 0.01;
    % num = ((finish - start) / increment) + 1;
    % C_T_array = zeros(1, num);
    % C_Q_array = zeros(1, num);
    % T_array = zeros(1, num);
    % Q_array = zeros(1, num);
    % J_array = start: increment: finish;

    % for i = 1: num
    %     J = J_array(1, i);
    %     U_1 = n * D * J;
    %     [C_T_array(1, i), C_Q_array(1, i), T_array(1, i), Q_array(1, i)] = vary(r_num, sigma_num, r_t, R, SIGMA, U_1, Omega, J, rho, B, D, n, theta, naca_matrix, sg_matrix, NACA6409, SG6043);
    %     fprintf('i=%d\n', i);
    % end
    % plot_graph(J_array, C_T_array, 'C_T-J', 'J', 'C_T');
    % plot_graph(J_array, C_Q_array, 'C_Q-J', 'J', 'C_Q');
    % plot_graph(J_array, T_array, 'T-J', 'J', 'C_Q');
    % plot_graph(J_array, Q_array, 'Q-J', 'J', 'C_Q');
    % file_name = 'parameters.xlsx';
    % writematrix(['J', 'C_T', 'C_Q'], file_name, 'Sheet', 2, 'Range', 'A1');
    % writematrix(([J_array; C_T_array; C_Q_array])', file_name, 'Sheet', 2, 'Range', 'A2');
    % writematrix([T Q], file_name, 'Sheet', 1, 'Range', 'B5')
    fprintf('Done\n');
end

function [C_T, C_Q, T, Q] = vary(r_num, sigma_num, r_t, R, SIGMA, U_1, Omega, J, rho, B, D, n, theta, naca_matrix, sg_matrix, NACA6409, SG6043)
    phi_matrix = NaN(r_num, sigma_num);
    a_matrix = NaN(r_num, sigma_num);
    b_matrix = NaN(r_num, sigma_num);
    C_L_matrix = NaN(r_num, sigma_num);
    C_D_matrix = NaN(r_num, sigma_num);
    for r = 1:r_num
        for sigma = 1:sigma_num
            if r < r_t
                [a, b, phi, C_L, C_D] = cal_a_b_phi(((r / r_num) * R), ((sigma / sigma_num) * SIGMA), U_1, R, Omega, J, SG6043, 1, sg_matrix, theta, r_num);
            else
                [a, b, phi, C_L, C_D] = cal_a_b_phi(((r / r_num) * R), ((sigma / sigma_num) * SIGMA), U_1, R, Omega, J, NACA6409, 1, naca_matrix, theta, r_num);
            end
            a_matrix(r, sigma) = a;
            b_matrix(r, sigma) = b;
            phi_matrix(r, sigma) = phi;
            C_L_matrix(r, sigma) = C_L;
            C_D_matrix(r, sigma) = C_D;
        end
    end
    r_arr = R / r_num:R / r_num:R;
    [phi_array, a_array, ~  , c_array, C_L_array, C_D_array] = pick_phi(phi_matrix, a_matrix, b_matrix, C_L_matrix, C_D_matrix, r_num, sigma_num, R, SIGMA);
    [C_T, C_Q, T, Q] = thrust_torque_coefficients(r_arr, rho, B, D, J, n, a_array, c_array, phi_array, C_L_array, C_D_array, r_num);
end



function [naca, sg] = get_C_L_C_D_matrix()
    naca = readtable('data\\naca6409.txt');
    sg = readtable('data\\sg6043.txt');
    naca = naca{:, :};
    sg = sg{:, :};
end


function [theta, naca_matrix, sg_matrix, NACA6409, SG6043, C_T, C_Q] = get_standard(r_num, sigma_num, r_t, R, U_1, Omega, J, rho, B, D, n, SIGMA)
    SG6043 = read_turbine_data('SG6043'); % 1 -> alpha, 2 -> C_L, 3 -> C_D
    NACA6409 = read_turbine_data('NACA6409');
    [naca_matrix, sg_matrix] = get_C_L_C_D_matrix();
    phi_matrix = NaN(r_num, sigma_num);
    a_matrix = NaN(r_num, sigma_num);
    b_matrix = NaN(r_num, sigma_num);

    for r = 1:r_num
        for sigma = 1:sigma_num
            if r < r_t
                [a, b, phi] = cal_a_b_phi(r / r_num * R, sigma / sigma_num * SIGMA, U_1, R, Omega, J, SG6043, 0);
            else
                [a, b, phi] = cal_a_b_phi(r / r_num * R, sigma / sigma_num * SIGMA, U_1, R, Omega, J, NACA6409, 0);
            end
            a_matrix(r, sigma) = a;
            b_matrix(r, sigma) = b;
            phi_matrix(r, sigma) = phi;
        end

    end

    r_arr = R / r_num:R / r_num:R;
    C_L_matrix = [(SG6043(2) + zeros(r_t, sigma_num)); (NACA6409(2) + zeros(r_num - r_t, sigma_num))];
    C_D_matrix = [(SG6043(3) + zeros(r_t, sigma_num)); (NACA6409(3) + zeros(r_num - r_t, sigma_num))];
    [phi_array, a_array, b_array, c_array, C_L_array, C_D_array] = pick_phi(phi_matrix, a_matrix, b_matrix, C_L_matrix, C_D_matrix, r_num, sigma_num, R, SIGMA);
    [C_T, C_Q, T, Q] = thrust_torque_coefficients(r_arr, rho, B, D, J, n, a_array, c_array, phi_array, C_L_array, C_D_array, r_num);
    theta = phi_array - [(zeros(1, r_t) + deg2rad(SG6043(1))), (zeros(1, r_num - r_t) + deg2rad(NACA6409(1)))];
    % plot_graph(r_arr, theta, 'theta -r', 'r (m)', 'theta (degree)');
end

function plot_graph(x, y, graph_title, x_axis, y_axis)
    figure
    plot(x, y, '-*')
    grid on
    title(graph_title)
    xlabel(x_axis)
    ylabel(y_axis)
end

function [C_T, C_Q, T, Q] = thrust_torque_coefficients(r, rho, B, D, J, n, a, c, phi, C_L, C_D, size)
    dT = NaN(1, size);
    dQ = NaN(1, size);
    M = C_L .* cos(phi);
    N = C_D .* sin(phi);
    something = rho * B * c .* ((1 - a) ./ sin(phi)).^2 * (J^2 * n^2 * D^2) / 2;
    dT = something .* (M + N);
    dQ = something .* r .* (M - N);
    Q = trapz(r, dQ);
    T = trapz(r, dT);
    something = B * rho * (n^2) * (D ^ 4);
    C_T = T / something;
    C_Q = Q / something / D;
end

function [phi_array, a_array, b_array, c_array, C_L_array, C_D_array] = pick_phi(phi_matrix, a_matrix, b_matrix, C_L_matrix, C_D_matrix, row, col, R, SIGMA)
    phi_array = NaN(1, row);
    a_array = NaN(1, row);
    b_array = NaN(1, row);
    c_array = NaN(1, row);
    C_L_array = NaN(1, row);
    C_D_array = NaN(1, row);
    i = 1;
    T = isnan(phi_matrix);
    for r = 1:row
        for c = 1:col
            if c < col - 1 && ~T(r, c + 1)
                continue;
            end
            phi_array(1, i) = phi_matrix(r, c);
            a_array(1, i) = a_matrix(r, c);
            b_array(1, i) = b_matrix(r, c);
            c_array(1, i) = (2 * pi * ((r / row) * R) * ((c / col) * SIGMA)) / 3  ;
            C_L_array(1, i) = C_L_matrix(r, c);
            C_D_array(1, i) = C_D_matrix(r, c);
            i = i + 1;
            break
        end
    end
end

function [B, U_1, R, D, Omega, rho, n, J, SIGMA] = read_data()
    text = fileread("data\common_data.txt");
    text = split(text);
    i = 1;

    while ~strcmp(text(i, 1), "-")

        if text(i, 1) == "B"
            B = str2double(text(i + 1, 1));
        elseif text(i, 1) == "U_1"
            U_1 = str2double(text(i + 1, 1));
        elseif text(i, 1) == "R"
            R = str2double(text(i + 1, 1));
        elseif text(i, 1) == "Omega"
            Omega = str2double(text(i + 1, 1));
        elseif text(i, 1) == "rho"
            rho = str2double(text(i + 1, 1));
        elseif text(i, 1) == "SIGMA"
            SIGMA = str2double(text(i + 1, 1));
        end

        i = i + 2;
    end

    D = 2 * R;
    n = Omega / 2 / pi;
    J = U_1 / n / D;
end

function turbine = read_turbine_data(turbine_name)
    turbine = NaN(3, 1);
    text = fileread("data\common_data.txt");
    text = split(text);
    i = 1;

    while ~strcmp(text(i, 1), turbine_name)
        i = i + 1;
    end

    i = i + 1;

    while 1

        if text(i, 1) == "alpha"
            turbine(1, 1) = str2double(text(i + 1, 1));
        elseif text(i, 1) == "C_L"
            turbine(2, 1) = str2double(text(i + 1, 1));
        elseif text(i, 1) == "C_D"
            turbine(3, 1) = str2double(text(i + 1, 1));
            break;
        end

        i = i + 2;
    end

end

function [a, b, phi, C_L, C_D] = cal_a_b_phi(r, sigma, U_1, R, Omega, J, turbine, vary, C_L_C_D_matrix, theta, r_num)
    if ~isfloat(r)
        error('r must be a double.')
    end

    if ~isfloat(sigma)
        error('sigma must be a double.')
    end

    if sigma < 0 || sigma > 1
        error('sigma must be in range [0, 0.5]. Input: %f', sigma)
    end

    D = 2 * R; a = 0.5; b = 0.5; i = 1; bound = 200;

    if r * sigma > 0.3
        bound = 10;
    end
    if vary == 0
        C_L = turbine(2, 1);
        C_D = turbine(3, 1);
    end
    while true
        pa = a; pb = b;
        phi = cal_phi(a, b, r, U_1, Omega);
        if vary == 1
            y = theta(1, floor((r / R) * r_num));
            [C_L, C_D] = cal_C_L_C_D(C_L_C_D_matrix, phi, y);
        end
        a = cal_a(phi, sigma, pa, C_L, C_D);
        b = cal_b(phi, sigma, pa, r, C_L, C_D, J, D);
        i = i + 1;
       
        
        if i > bound
            a = NaN; b = NaN; phi = NaN; C_L = NaN; C_D = NaN;
            break;
        elseif abs(a - pa) + abs(b - pb) <= 0.0001
            break;
        end

    end

end

function [C_L, C_D] = cal_C_L_C_D(matrix, phi, theta)
    alpha = (phi - theta) * 180 / pi;
    index = floor((alpha - matrix(1, 1)) * 4);
    ratio = ((alpha - matrix(1, 1)) * 4) - index;

    if 1 <= index && index <= 95
        C_L = (1 - ratio) * matrix(index, 2) + ratio * matrix(index + 1, 2);
        C_D = (1 - ratio) * matrix(index, 3) + ratio * matrix(index + 1, 3);
    elseif index == 96
        C_L = (1 - ratio) * matrix(index, 2);
        C_D = (1 - ratio) * matrix(index, 3);
    elseif index > 96
        C_L = matrix(96, 2);
        C_D = matrix(96, 3);
    elseif index < 1
        C_L = matrix(1, 2);
        C_D = matrix(1, 3);
        %C_L = NaN;
        %C_D = NaN;
    end

end

function phi = cal_phi(a, b, r, U_1, Omega)
    phi = atan(U_1 * (1 - a) / (Omega * r * (1 + b)));
end

function a = cal_a(phi, sigma, pa, C_L, C_D)
    a = sigma * (1 - pa) * (C_L * cos(phi) + C_D * sin(phi)) / (4 * (sin(phi))^2);
end

function b = cal_b(phi, sigma, pa, r, C_L, C_D, J, D)
    b = sigma * (1 - pa) * J * (C_L * sin(phi) - C_D * cos(phi)) / (4 * (sin(phi))^2 * 2 * pi * r / D);
end
