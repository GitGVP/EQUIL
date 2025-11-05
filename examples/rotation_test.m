% Actually testing all the git stuff too
% Testing the rotation

[Le2, LXe2] = equilSol('beta',0.5, 's0', 5, 'debug', 4,...
    'residuals_fun',@residuals_rotation,'jacobian_fun', @jacobian_rotation, ...
    'equation_of_state', @isotropic,'Ns', 2);
LXe2.eps_val = 0.37225;


LXe2.Sbc(1) = 0;
LXe2.S2bc = 0;
LXe2.Sbc(2) = 0;
LXe2.S3bc = 0;
LYe2 = equilY(Le2, LXe2);