%% --- CORAL function (same as before) ---
function Xs_new = coral(Xs, Xt)
    Cs = cov(Xs) + eye(size(Xs,2))*1e-6;
    Ct = cov(Xt) + eye(size(Xt,2))*1e-6;

    [Us, Ss] = eig(Cs);
    [Ut, St] = eig(Ct);

    Ds_inv_sqrt = diag(diag(Ss).^(-0.5));
    Dt_sqrt     = diag(diag(St).^(+0.5));

    whitening = Us * Ds_inv_sqrt * Us';
    coloring  = Ut * Dt_sqrt   * Ut';

    Xs_white = (whitening * Xs')';
    Xs_new   = (coloring  * Xs_white')';

    % % Align means
    % mu_s  = mean(Xs_new,1);
    % mu_t  = mean(Xt,1);
    % Xs_new = bsxfun(@minus, Xs_new, mu_s) + mu_t;
end
