% exact solution for 2D model

function exact_2d(param_W, init_s, fID)
    tic;

    N = 2;
    L = 28; % 32;
    M = 256; % 384; 

    param_C = 0.25;
    %param_W = 10.0;

    % angle to rotate
    rot = @(r,theta) [r(1)*cos(theta) + r(2)*sin(theta), -r(1)*sin(theta) + r(2)*cos(theta)];
    theta = 0.5 * pi;

    theta1 = 0.0;
    theta2 = theta;

    % initial position, momentum, and diab surface
    r0 = [-9, 0];
    k0 = [18, 0];
    r01 = rot(r0,-theta1);
    k01 = rot(k0,-theta1);
    r02 = rot(r0,-theta2);
    k02 = rot(k0,-theta2);

    if init_s == 0
        r = r01;
        k = k01;
        init_s = 0;
    elseif init_s == 1
        r = r02;
        k = k02;
        init_s = 1;
    else
        fslkfjlkas
    end

    % setup
    xI = r(1);
    yI = r(2);
    kxI = k(1);
    kyI = k(2);
    sigmax = 1;
    sigmay = 1;

    c1 = sqrt(1 - init_s);
    c2 = sqrt(init_s);

    dt = 0.2;
    Nstep = 100000; 
    output_step = 2500;
    mass = 1000.0;

    enable_plot = false;
    enable_abc = true;

    % potential energy 
    debug_flat_potential = false;
    if debug_flat_potential == true
        fprintf(fID, '# flat potential, for debug\n');
        debug_flat_potential = true;
        cal_H11 = @(r) -0.02 * cos(0.5 * pi * (erf(3 * r(1)) + 1));
        cal_H22 = @(r) 0.02 * cos(0.5 * pi * (erf(3 * r(1)) + 1));
        cal_H12 = @(r) 0.02 * sin(0.5 * pi * (erf(3 * r(1)) + 1)) * exp(1i * param_W * r(2));
    else
        fprintf(fID, '# conner potential\n');
        cal_H_base = @(r) tanh(r(2)-5) - tanh(r(2)+5) + tanh(r(1)) + 3;
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

    % for abc
    if enable_abc == true
        U0 = 0.04;
        alpha = 0.1;
        reduce_xL = zeros(M,M);
        reduce_xR = zeros(M,M);
        reduce_yU = zeros(M,M);
        reduce_yD = zeros(M,M);
        for j=1:M
            for k=1:M
                reduce_xL(k,j) = (1.0 - U0 / cosh(alpha * (j-1))^2 * dt);
                reduce_xR(k,j) = (1.0 - U0 / cosh(alpha * (M-j+1))^2 * dt);
                reduce_yD(k,j) = (1.0 - U0 / cosh(alpha * (k-1))^2 * dt);
                reduce_yU(k,j) = (1.0 - U0 / cosh(alpha * (M-k+1))^2 * dt);
            end
        end
        xL_absorb_accu1 = 0.0;
        xR_absorb_accu1 = 0.0;
        xL_absorb_accu2 = 0.0;
        xR_absorb_accu2 = 0.0;
        yD_absorb_accu1 = 0.0;
        yD_absorb_accu2 = 0.0;
        yU_absorb_accu1 = 0.0;
        yU_absorb_accu2 = 0.0;
    end
    % construct TU on k grid
    T = (meshkx.^2 + meshky.^2) / 2 / mass;
    TU = exp(-1i * dt * T);
    TU = fftshift(TU);
    % construct VU
    VU = zeros(M,M,N,N);
    Hs = zeros(M,M,N,N);
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
            % VU & Hs
            VU(k,j,:,:) = expm(-1i * dt / 2 * H);
            Hs(k,j,:,:) = H;
        end
    end
    % Initial wavefunction -- Gaussian wavepacket 
    psi0 = zeros(M,M,N);
    psi0(:,:,1) = c1 * exp(1i*(kxI*meshx + kyI*meshy)) .* exp(-(meshx-xI).^2/sigmax^2 - (meshy-yI).^2/sigmay^2);
    psi0(:,:,2) = c2 * exp(1i*(kxI*meshx + kyI*meshy)) .* exp(-(meshx-xI).^2/sigmax^2 - (meshy-yI).^2/sigmay^2);
    psi0 = psi0 / sqrt(sum(sum(sum(abs(psi0).^2))));
    % psim & psi_k_m -- for plot
    psi0_k(:,:,1) = fftshift(fft2(psi0(:,:,1)));
    psi0_k(:,:,2) = fftshift(fft2(psi0(:,:,2)));
    psi0_k = psi0_k / sqrt(sum(sum(sum(abs(psi0_k).^2))));
    psi_max = max(max(max(abs(psi0).^2)));
    psi_k_max = max(max(max(abs(psi0_k).^2)));
    % propagate WF
    ana_step = 0;
    psi = psi0;
    for t=0:Nstep
        % exp(-iVdt/2) * |Psi>
        psi_k(:,:,1) = VU(:,:,1,1).*psi(:,:,1) + VU(:,:,1,2).*psi(:,:,2);
        psi_k(:,:,2) = VU(:,:,2,1).*psi(:,:,1) + VU(:,:,2,2).*psi(:,:,2);
        % exp(-iTdt) * psi
        psi_k(:,:,1) = TU .* fft2(psi_k(:,:,1));
        psi_k(:,:,2) = TU .* fft2(psi_k(:,:,2));
        % exp(-iVdt/2) * psi
        psi_k(:,:,1) = ifft2(psi_k(:,:,1));
        psi_k(:,:,2) = ifft2(psi_k(:,:,2));
        psi(:,:,1) = VU(:,:,1,1).*psi_k(:,:,1) + VU(:,:,1,2).*psi_k(:,:,2);
        psi(:,:,2) = VU(:,:,2,1).*psi_k(:,:,1) + VU(:,:,2,2).*psi_k(:,:,2);
        % abosrbing boundary condition, notice energy conservation will break if abc is enabled
        if enable_abc == true
            % reduce xL boundary 
            n1_before_reduce = sum(sum(abs(psi(:,:,1)).^2));
            n2_before_reduce = sum(sum(abs(psi(:,:,2)).^2));
            psi(:,:,1) = psi(:,:,1) .* reduce_xL;
            psi(:,:,2) = psi(:,:,2) .* reduce_xL;
            n1_after_reduce = sum(sum(abs(psi(:,:,1)).^2));
            n2_after_reduce = sum(sum(abs(psi(:,:,2)).^2));
            xL_absorb_accu1 = xL_absorb_accu1 + n1_before_reduce - n1_after_reduce;
            xL_absorb_accu2 = xL_absorb_accu2 + n2_before_reduce - n2_after_reduce;
            % reduce xR boundary 
            n1_before_reduce = sum(sum(abs(psi(:,:,1)).^2));
            n2_before_reduce = sum(sum(abs(psi(:,:,2)).^2));
            psi(:,:,1) = psi(:,:,1) .* reduce_xR;
            psi(:,:,2) = psi(:,:,2) .* reduce_xR;
            n1_after_reduce = sum(sum(abs(psi(:,:,1)).^2));
            n2_after_reduce = sum(sum(abs(psi(:,:,2)).^2));
            xR_absorb_accu1 = xR_absorb_accu1 + n1_before_reduce - n1_after_reduce;
            xR_absorb_accu2 = xR_absorb_accu2 + n2_before_reduce - n2_after_reduce;
            % reduce yD boundary
            n1_before_reduce = n1_after_reduce;
            n2_before_reduce = n2_after_reduce;
            psi(:,:,1) = psi(:,:,1) .* reduce_yD;
            psi(:,:,2) = psi(:,:,2) .* reduce_yD;
            n1_after_reduce = sum(sum(abs(psi(:,:,1)).^2));
            n2_after_reduce = sum(sum(abs(psi(:,:,2)).^2));
            yD_absorb_accu1 = yD_absorb_accu1 + n1_before_reduce - n1_after_reduce;
            yD_absorb_accu2 = yD_absorb_accu2 + n2_before_reduce - n2_after_reduce;
            % reduce yU boundary
            n1_before_reduce = n1_after_reduce;
            n2_before_reduce = n2_after_reduce;
            psi(:,:,1) = psi(:,:,1) .* reduce_yU;
            psi(:,:,2) = psi(:,:,2) .* reduce_yU;
            n1_after_reduce = sum(sum(abs(psi(:,:,1)).^2));
            n2_after_reduce = sum(sum(abs(psi(:,:,2)).^2));
            yU_absorb_accu1 = yU_absorb_accu1 + n1_before_reduce - n1_after_reduce;
            yU_absorb_accu2 = yU_absorb_accu2 + n2_before_reduce - n2_after_reduce;
        end
        % analysis & report
        if mod(t,output_step) == 0
            % header
            if t == 0
                fprintf(fID, '# EXACT DIAB\n');
                fprintf(fID, '# xI = %8.4f yI = %8.4f kxI = %8.4f kyI = %8.4f sigmax = %8.4f sigmay = %8.4f W = %8.4f init_s = %8.4f c1 = %8.4f c2 = %8.4f \n', ...
                                    xI, yI, kxI, kyI, sigmax, sigmay, param_W, init_s, c1, c2);
                fprintf(fID, '# L = %8.4f M = %8d dt = %8.4f Nstep = %8d output_step = %8d\n', ...
                                    L, M, dt, Nstep, output_step);
                if enable_abc == true
                    fprintf(fID, '#%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s\n', ...
                                    't', 'n1d', 'n2d', 'Etot', ...
                                    'xL_absorb1', 'xL_absorb2', 'xR_absorb1', 'xR_absorb2',  ...
                                    'yD_absorb1', 'yD_aborb2', 'yU_absorb1', 'yU_aborb2', ...
                                    'tot_absorb1', 'tot_absorb2', 'tot_absorb');
                else
                    fprintf(fID, '#%16s%16s%16s%16s%16s%16s\n', ...
                                    't', 'n1d', 'n2d', 'KE', 'PE', 'Etot');
                end
            end
            % analysis
            psi_k(:,:,1) = fftshift(fft2(psi(:,:,1)));
            psi_k(:,:,2) = fftshift(fft2(psi(:,:,2)));
            psi_k = psi_k / sqrt(sum(sum(sum(abs(psi_k).^2))));
            nd1 = sum(sum(abs(psi(:,:,1)).^2));
            nd2 = sum(sum(abs(psi(:,:,2)).^2));
            KE = sum(sum((abs(psi_k(:,:,1)).^2 + abs(psi_k(:,:,2)).^2) .* (meshkx.^2 + meshky.^2) / 2 / mass));
            PE = sum(sum( conj(psi(:,:,1)) .* Hs(:,:,1,1) .* psi(:,:,1) + conj(psi(:,:,1)) .* Hs(:,:,1,2) .* psi(:,:,2) + conj(psi(:,:,2)) .* Hs(:,:,2,1) .* psi(:,:,1) + conj(psi(:,:,2)) .* Hs(:,:,2,2) .* psi(:,:,2) ));
            if enable_abc == true
                fprintf(fID, '#%16.10f%16.10f%16.10f%16.10f%16.10f%16.10f%16.10f%16.10f%16.10f%16.10f%16.10f%16.10f%16.10f%16.10f%16.10f\n', ...
                            t*dt, ...
                            nd1, nd2, ...
                            KE+PE, ...
                            xL_absorb_accu1, xL_absorb_accu2, ...
                            xR_absorb_accu1, xR_absorb_accu2, ...
                            yD_absorb_accu1, yD_absorb_accu2, ...
                            yU_absorb_accu1, yU_absorb_accu2, ...
                            xL_absorb_accu1+xR_absorb_accu1+yD_absorb_accu1+yU_absorb_accu1, ...
                            xL_absorb_accu2+xR_absorb_accu2+yD_absorb_accu2+yU_absorb_accu2, ...
                            xL_absorb_accu1+xR_absorb_accu1+yD_absorb_accu1+yU_absorb_accu1+xL_absorb_accu2+xR_absorb_accu2+yD_absorb_accu2+yU_absorb_accu2 ...
                            );
            else
                fprintf(fID, '#%16.10f%16.10f%16.10f%16.10f%16.10f%16.10f\n', ...
                            t*dt, ...
                            nd1, nd2, ...
                            KE, PE, KE+PE ...
                            );
            end
            % plot
            if enable_plot == true
                subplot(2,2,1);
                contour(y0,x0,abs(psi(:,:,1)).^2,[psi_max/100:psi_max/100:psi_max]);
                title('Real space -- Pop Diab 1');

                subplot(2,2,2);
                contour(y0,x0,abs(psi(:,:,2)).^2,[psi_max/100:psi_max/100:psi_max]);
                title('Real space -- Pop Diab 2');

                subplot(2,2,3);
                contour(ky0,kx0,abs(psi_k(:,:,1)).^2,[psi_k_max/100:psi_k_max/100:psi_k_max]);
                title('Mom space -- Pop Diab 1');

                subplot(2,2,4);
                contour(ky0,kx0,abs(psi_k(:,:,2)).^2,[psi_k_max/100:psi_k_max/100:psi_k_max]); 
                title('Mom space -- Pop Diab 2');

                drawnow;
            end
            % check ending
            if enable_abc == true
                threash = 0.999;
                if xL_absorb_accu1+xR_absorb_accu1+yD_absorb_accu1+yU_absorb_accu1+xL_absorb_accu2+xR_absorb_accu2+yD_absorb_accu2+yU_absorb_accu2 > threash
                    fprintf('# %.3f pop has been absorbed, ending\n', threash);
                    break;
                end
            end
        end
    end
    fprintf('# '); toc; fprintf('\n');
    % final output
    if enable_abc == true
        fprintf(fID, ' %16.10f%16.10f\n', ...
                    xL_absorb_accu1+xR_absorb_accu1+yD_absorb_accu1+yU_absorb_accu1, ...
                    xL_absorb_accu2+xR_absorb_accu2+yD_absorb_accu2+yU_absorb_accu2 ...
                    );
    end
end
