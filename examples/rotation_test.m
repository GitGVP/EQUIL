% Actually testing all the git stuff too
% Testing the rotation

[L, LX] = equilSol('beta',0.5, 's0', 5, 'debug', 4,...
    'residuals_fun',@residuals_rotation,'jacobian_fun', @jacobian_rotation, ...
    'equation_of_state', @parallel_rotating,'Ns', 2);
LX.eps_val = 0.37225;


LX.Sbc(1) = 0;
LX.S2bc = 0;
LX.Sbc(2) = 0;
LX.S3bc = 0;
LY = equilY(L, LX);