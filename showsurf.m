
function showsurf()
    L = 24.0;
    M = 64;

    debug_flat_potential = false;

    if debug_flat_potential == true
        fprintf('# flat potential, for debug\n');
        debug_flat_potential = true;
        param_W = 0.5;
        cal_H11 = @(x,y) -0.02 * cos(0.5 * pi * (erf(3 * x) + 1));
        cal_H22 = @(x,y) 0.02 * cos(0.5 * pi * (erf(3 * x) + 1));
        cal_H12 = @(x,y) 0.02 * sin(0.5 * pi * (erf(3 * x) + 1)) * exp(1i * param_W * y);
    else
        fprintf('# conner potential\n');
        param_C = 0.1;
        param_W = 0.5;
        cal_H11 = @(x,y) tanh(x-5) - tanh(x+5) + tanh(y) + 3; % + 0.1 * 5^2 / ((x-y)^2 + 5^2);
        cal_H22 = @(x,y) tanh(y-5) - tanh(y+5) + tanh(x) + 3; % + 0.1 * 5^2 / ((x-y)^2 + 5^2);
        cal_H12 = @(x,y) param_C * exp(1i * param_W * sqrt(x^2 + y^2));
    end

    % grids
    x0 = linspace(-L/2, L/2, M)';
    y0 = linspace(-L/2, L/2, M)';
    dx = x0(2) - x0(1);
    dy = y0(2) - y0(1);
    dkx = 2 * pi / M / dx;
    dky = 2 * pi / M / dy;
    kx0 = (-M/2:M/2-1)' * dkx;
    ky0 = (-M/2:M/2-1)' * dky;
    [meshx, meshy] = meshgrid(x0, y0);
    [meshkx, meshky] = meshgrid(kx0, ky0);


    H11 = zeros(M,M);
    H22 = zeros(M,M);
    E1 = zeros(M,M);
    E2 = zeros(M,M);
    for j=1:M
        for k=1:M
            x = x0(j);
            y = y0(k);

            H = zeros(2,2);
            H(1,1) = cal_H11(x,y);
            H(2,2) = cal_H22(x,y);
            H(1,2) = cal_H12(x,y);
            H(2,1) = conj(H(1,2));

            eva = eig(H);
            H11(j,k) = H(1,1);
            H22(j,k) = H(2,2);
            E1(j,k) = eva(1);
            E2(j,k) = eva(2);
        end
    end

    figure;
    subplot(2,2,1);
    surf(meshx, meshy, H11);
    xlabel('x');
    ylabel('y');
    title('H11');

    subplot(2,2,2);
    surf(meshx, meshy, H22);
    xlabel('x');
    ylabel('y');
    title('H22');

    subplot(2,2,3);
    surf(meshx, meshy, E1);
    xlabel('x');
    ylabel('y');
    title('E1');

    subplot(2,2,4);
    surf(meshx, meshy, E2);
    xlabel('x');
    ylabel('y');
    title('E2');

    drawnow;
end
