% exact solution for 2D model

function exact_2d(r, k, init_s, fID)
    tic;

    xI = r(1);
    yI = r(2);
    kxI = k(1);
    kyI = k(2);
    sigmax = 1;
    sigmay = 1;

    c1 = sqrt(1 - init_s);
    c2 = sqrt(init_s);

    N = 2;
    L = 24; % 32;
    M = 192; % 384; 

    dt = 0.2;
    Nstep = 50000; 
    output_step = 500;
    mass = 1000.0;

    enable_plot = false;
    enable_abc = true;

    % potential energy 
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
        param_C = 0.15;
        param_W = 0.0;
        %{
        cal_H11 = @(x,y) tanh(x-5) - tanh(x+5) + tanh(y) + 3; % + 0.1 * 5^2 / ((x-y)^2 + 5^2);
        cal_H22 = @(x,y) tanh(y-5) - tanh(y+5) + tanh(x) + 3; % + 0.1 * 5^2 / ((x-y)^2 + 5^2);
        cal_H12 = @(x,y) param_C * exp(1i * param_W * sqrt(x^2 + y^2));
        %}
        cal_H11 = @(x,y) tanh(x-8) - tanh(x+2) + tanh(y-3) + 3; % + 0.1 * 5^2 / ((x-y)^2 + 5^2);
        cal_H22 = @(x,y) tanh(y-8) - tanh(y+2) + tanh(x-3) + 3; % + 0.1 * 5^2 / ((x-y)^2 + 5^2);
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
    % for abc
    if enable_abc == true
        U0 = 0.04;
        alpha = 0.1;
        reduce_x = zeros(M,M);
        reduce_y = zeros(M,M);
        for j=1:M
            for k=1:M
                reduce_x(k,j) = (1.0 - U0 / cosh(alpha * (j-1))^2 * dt) ...
                                    * (1.0 - U0 / cosh(alpha * (M-j+1))^2 * dt);
                reduce_y(k,j) = (1.0 - U0 / cosh(alpha * (k-1))^2 * dt) ...
                                     * (1.0 - U0 / cosh(alpha * (M-k+1))^2 * dt);
            end
        end
        x_absorb_accu1 = 0.0;
        x_absorb_accu2 = 0.0;
        y_absorb_accu1 = 0.0;
        y_absorb_accu2 = 0.0;
    end
    % construct TU on k grid
    T = (meshkx.^2 + meshky.^2) / 2 / mass;
    TU = exp(-1i * dt * T);
    TU = fftshift(TU);
    % construct VU
    VU = zeros(M,M,N,N);
    Hs = zeros(M,M,N,N);
    evas = zeros(M,M,N,N);
    evts = zeros(M,M,N,N);
    for j=1:M
        for k=1:M
            x = x0(j);
            y = y0(k);
            % construct H
            H = zeros(2,2);
            H(1,1) = cal_H11(x,y);
            H(2,2) = cal_H22(x,y);
            H(1,2) = cal_H12(x,y);
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
            % VU & Hs
            VU(k,j,:,:) = expm(-1i * dt / 2 * H);
            Hs(k,j,:,:) = H;
        end
    end
    % Initial wavefunction -- Gaussian wavepacket on adiab
    psiad0 = zeros(M,M,N);
    psiad0(:,:,1) = c1 * exp(1i*(kxI*meshx + kyI*meshy)) .* exp(-(meshx-xI).^2/sigmax^2 - (meshy-yI).^2/sigmay^2);
    psiad0(:,:,2) = c2 * exp(1i*(kxI*meshx + kyI*meshy)) .* exp(-(meshx-xI).^2/sigmax^2 - (meshy-yI).^2/sigmay^2);
    psiad0 = psiad0 / sqrt(sum(sum(sum(abs(psiad0).^2))));
    % convert to diab for propagation
    psi0 = zeros(M,M,N);
    psi0(:,:,1) = evts(:,:,1,1) .*  psiad0(:,:,1) + evts(:,:,1,2) .* psiad0(:,:,2);
    psi0(:,:,2) = evts(:,:,2,1) .*  psiad0(:,:,1) + evts(:,:,2,2) .* psiad0(:,:,2);
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
    for t=0:Nstep-1
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
            % reduce x boundary 
            n1_before_reduce = sum(sum(abs(psi(:,:,1)).^2));
            n2_before_reduce = sum(sum(abs(psi(:,:,2)).^2));

            psi(:,:,1) = psi(:,:,1) .* reduce_x;
            psi(:,:,2) = psi(:,:,2) .* reduce_x;

            n1_after_reduce = sum(sum(abs(psi(:,:,1)).^2));
            n2_after_reduce = sum(sum(abs(psi(:,:,2)).^2));

            x_absorb_accu1 = x_absorb_accu1 + n1_before_reduce - n1_after_reduce;
            x_absorb_accu2 = x_absorb_accu2 + n2_before_reduce - n2_after_reduce;

            % reduce y boundary
            n1_before_reduce = n1_after_reduce;
            n2_before_reduce = n2_after_reduce;

            psi(:,:,1) = psi(:,:,1) .* reduce_y;
            psi(:,:,2) = psi(:,:,2) .* reduce_y;

            n1_after_reduce = sum(sum(abs(psi(:,:,1)).^2));
            n2_after_reduce = sum(sum(abs(psi(:,:,2)).^2));

            y_absorb_accu1 = y_absorb_accu1 + n1_before_reduce - n1_after_reduce;
            y_absorb_accu2 = y_absorb_accu2 + n2_before_reduce - n2_after_reduce;
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
                    fprintf(fID, '#%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s\n', ...
                                    't', 'n1d', 'n2d', 'KE', 'PE', 'Etot', ...
                                    'x_absorb1', 'x_absorb2', 'y_absorb1', 'y_aborb2', 'tot_absorb');
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
                fprintf(fID, '#%16.10f%16.10f%16.10f%16.10f%16.10f%16.10f%16.10f%16.10f%16.10f%16.10f%16.10f\n', ...
                            t*dt, ...
                            nd1, nd2, ...
                            KE, PE, KE+PE, ...
                            x_absorb_accu1, x_absorb_accu2, ...
                            y_absorb_accu1, y_absorb_accu2, ...
                            x_absorb_accu1+x_absorb_accu2+y_absorb_accu1+y_absorb_accu2...
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
                if x_absorb_accu1+x_absorb_accu2+y_absorb_accu1+y_absorb_accu2 > threash
                    fprintf('# %.3f pop has been absorbed, ending\n', threash);
                    break;
                end
            end
        end
    end
    fprintf('# '); toc; fprintf('\n');
end
