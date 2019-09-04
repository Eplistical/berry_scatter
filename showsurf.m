
function showsurf()
    N = 2;
    L = 30;
    M = 64;

    rot = @(r,theta) [r(1)*cos(theta) + r(2)*sin(theta), -r(1)*sin(theta) + r(2)*cos(theta)];
    theta =  0.5 * pi;

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
        param_C = 0.25;
        param_W = 0.0;

        theta1 = 0;
        theta2 = theta;

        r0 = [-10, 0];
        k0 = [18, 0];
        r01 = rot(r0,-theta1);
        k01 = rot(k0,-theta1);
        r02 = rot(r0,-theta2);
        k02 = rot(k0,-theta2);

        cal_H_base = @(r) tanh(r(2)-5) - tanh(r(2)+5) + tanh(r(1)) + 3; % + 0.1 * 5^2 / ((x-y)^2 + 5^2);
        cal_H11 = @(r) cal_H_base(rot(r,theta1));
        cal_H22 = @(r) cal_H_base(rot(r,theta2));
        cal_H12 = @(r) param_C * exp(1i * param_W * sqrt(r(1)^2 + r(2)^2));
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
    % ==> meshx(k,j) = x(j), meshy(k,j) = y(k)

    H11 = zeros(M,M);
    H22 = zeros(M,M);
    E1 = zeros(M,M);
    E2 = zeros(M,M);
    evas = zeros(M,M,N,N);
    evts = zeros(M,M,N,N);
    for j=1:M
        for k=1:M
            x = meshx(k,j);
            y = meshy(k,j);
            r = [x,y];
            % construct H
            H = zeros(2,2);
            H(1,1) = cal_H11(r);
            H(2,2) = cal_H22(r);
            H(1,2) = cal_H12(r);
            H(2,1) = conj(H(1,2));
            % eva & evt
            [evt, eva] = eig(H);
            if j == 1 && k == 1
                phase1 = evt(1,1) / abs(evt(1,1));
                phase2 = evt(2,2) / abs(evt(2,2));
                evt(:,1) = evt(:,1) * conj(phase1);
                evt(:,2) = evt(:,2) * conj(phase2);
            elseif k == 1
                phase1 = conj(evts(k,j-1,1,1)) * evt(1,1) + conj(evts(k,j-1,2,1)) * evt(2,1);
                phase2 = conj(evts(k,j-1,1,2)) * evt(1,2) + conj(evts(k,j-1,2,2)) * evt(2,2);
                evt(:,1) = evt(:,1) * conj(phase1);
                evt(:,2) = evt(:,2) * conj(phase2);
            else
                phase1 = conj(evts(k-1,j,1,1)) * evt(1,1) + conj(evts(k-1,j,2,1)) * evt(2,1);
                phase2 = conj(evts(k-1,j,1,2)) * evt(1,2) + conj(evts(k-1,j,2,2)) * evt(2,2);
                evt(:,1) = evt(:,1) * conj(phase1);
                evt(:,2) = evt(:,2) * conj(phase2);
            end
            evts(k,j,:,:) = evt;
            evas(k,j,:,:) = eva;
            % for plot 
            H11(k,j) = H(1,1);
            H22(k,j) = H(2,2);
            E1(k,j) = eva(1,1);
            E2(k,j) = eva(2,2);
        end
    end

    Eonly = false;
    enable3d = false;
    if Eonly == true
        lv = 10;

        figure;
        contourf(meshx', meshy', E2', lv);
        xlabel('x');
        ylabel('y');
        title('E2');
        hold on;
        scatter(r01(1), r01(2), 'g');
        plot([r01(1), r01(1) + k01(1) * 0.2], ...
             [r01(2), r01(2) + k01(2) * 0.2], 'g');
        scatter(r02(1), r02(2), 'r');
        plot([r02(1), r02(1) + k02(1) * 0.2], ...
             [r02(2), r02(2) + k02(2) * 0.2], 'r');

        figure;
        contourf(meshx, meshy, E1, lv);
        xlabel('x');
        ylabel('y');
        title('E1');
        hold on;
        scatter(r01(1), r01(2), 'g');
        plot([r01(1), r01(1) + k01(1) * 0.2], ...
             [r01(2), r01(2) + k01(2) * 0.2], 'g');
        scatter(r02(1), r02(2), 'r');
        plot([r02(1), r02(1) + k02(1) * 0.2], ...
             [r02(2), r02(2) + k02(2) * 0.2], 'r');
    else 
        if enable3d == true
            figure();

            subplot(2,2,1);
            surf(meshx', meshy', H11');
            xlabel('x');
            ylabel('y');
            title('H11');

            subplot(2,2,2);
            surf(meshx', meshy', H22');
            xlabel('x');
            ylabel('y');
            title('H22');

            subplot(2,2,3);
            surf(meshx', meshy', E1');
            xlabel('x');
            ylabel('y');
            title('E1');

            subplot(2,2,4);
            surf(meshx', meshy', E2');
            xlabel('x');
            ylabel('y');
            title('E2');
        else
            figure;
            lv = 10;

            subplot(2,2,1);
            contourf(meshx', meshy', H11', lv);
            xlabel('x');
            ylabel('y');
            title('H11');

            subplot(2,2,2);
            contourf(meshx', meshy', H22', lv);
            xlabel('x');
            ylabel('y');
            title('H22');

            subplot(2,2,3);
            contourf(meshx', meshy', E1', lv);
            xlabel('x');
            ylabel('y');
            title('E1');

            subplot(2,2,4);
            contourf(meshx', meshy', E2', lv);
            xlabel('x');
            ylabel('y');
            title('E2');
        end
    end


    drawnow;
end
