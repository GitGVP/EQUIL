function [Nval, N1, N2, N3] = bspline_eval_all(T, p, x)
    x = x(:)';                     % row vector
    nx = numel(x);
    K = length(T);
    nb = K - p - 1;

    % N(i,:,d+1) holds degree d basis i (shape nb x nx x (p+1))
    N = zeros(nb, nx, p+1);

    % degree 0
    for i = 1:nb
        N(i, :, 1) = (x >= T(i) & x < T(i+1));
    end
    % include right endpoint
    idxR = find(x == T(end));
    if ~isempty(idxR)
        N(:, idxR, 1) = 0;
        N(nb, idxR, 1) = 1;
    end

    % Cox-de Boor: build up to degree p
    for d = 1:p
        for i = 1:nb
            a = T(i+d)   - T(i);
            b = T(i+d+1) - T(i+1);
            term1 = zeros(1,nx);
            term2 = zeros(1,nx);
            if a > 0
                term1 = ((x - T(i)) / a) .* N(i, :, d);
            end
            if b > 0 && (i+1) <= nb
                term2 = ((T(i+d+1) - x) / b) .* N(i+1, :, d);
            end
            N(i, :, d+1) = term1 + term2;
        end
    end

    % degree-p values
    Nval = squeeze(N(:, :, p+1));   % nb x nx

    % first derivative N'_{i,p} = p/(t_{i+p}-t_i) N_{i,p-1} - p/(t_{i+p+1}-t_{i+1}) N_{i+1,p-1}
    N1 = zeros(nb, nx);
    for i = 1:nb
        a = T(i+p)   - T(i);
        b = T(i+p+1) - T(i+1);
        term1 = zeros(1,nx);
        term2 = zeros(1,nx);
        if a > 0
            term1 = (p / a) * N(i, :, p);
        end
        if b > 0 && (i+1) <= nb
            term2 = (p / b) * N(i+1, :, p);
        end
        N1(i, :) = term1 - term2;
    end

    % second derivative via N'_{i,p-1}
    if p >= 2
        Np1 = zeros(nb, nx); % N' for degree p-1
        for i = 1:nb
            a2 = T(i+p-1) - T(i);
            b2 = T(i+p)   - T(i+1);
            u1 = zeros(1,nx);
            u2 = zeros(1,nx);
            if a2 > 0
                u1 = ((p-1) / a2) * N(i, :, p-1);
            end
            if b2 > 0 && (i+1) <= nb
                u2 = ((p-1) / b2) * N(i+1, :, p-1);
            end
            Np1(i, :) = u1 - u2;
        end
        N2 = zeros(nb, nx);
        for i = 1:nb
            a = T(i+p)   - T(i);
            b = T(i+p+1) - T(i+1);
            t1 = zeros(1,nx);
            t2 = zeros(1,nx);
            if a > 0
                t1 = (p / a) * Np1(i, :);
            end
            if b > 0 && (i+1) <= nb
                t2 = (p / b) * Np1(i+1, :);
            end
            N2(i, :) = t1 - t2;
        end
    else
        N2 = zeros(nb, nx);
    end
    if p >= 3
        % First derivative of degree p-2
        D1 = zeros(nb, nx);
        for i = 1:nb
            a = T(i+p-2) - T(i);
            b = T(i+p-1) - T(i+1);
            t1 = zeros(1,nx);
            t2 = zeros(1,nx);
            if a > 0
                t1 = ((p-2) / a) * N(i, :, p-2);
            end
            if b > 0 && (i+1) <= nb
                t2 = ((p-2) / b) * N(i+1, :, p-2);
            end
            D1(i, :) = t1 - t2;
        end
    
        % Second derivative of degree p-1
        D2 = zeros(nb, nx);
        for i = 1:nb
            a = T(i+p-1) - T(i);
            b = T(i+p)   - T(i+1);
            t1 = zeros(1,nx);
            t2 = zeros(1,nx);
            if a > 0
                t1 = ((p-1) / a) * D1(i, :);
            end
            if b > 0 && (i+1) <= nb
                t2 = ((p-1) / b) * D1(i+1, :);
            end
            D2(i, :) = t1 - t2;
        end
    
        % Third derivative of degree p
        N3 = zeros(nb, nx);
        for i = 1:nb
            a = T(i+p)   - T(i);
            b = T(i+p+1) - T(i+1);
            t1 = zeros(1,nx);
            t2 = zeros(1,nx);
            if a > 0
                t1 = (p / a) * D2(i, :);
            end
            if b > 0 && (i+1) <= nb
                t2 = (p / b) * D2(i+1, :);
            end
            N3(i, :) = t1 - t2;
        end
    else
        N3 = zeros(nb, nx);
    end
end
