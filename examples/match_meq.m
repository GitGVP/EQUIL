%% Reproduce an analytic MEQ equilibrium with variational EQUIL

% add the path to MEQ and run with iterq >0 (q is necessary) and a higher
% resolution greatly improves the visual match.

meqpath = '/home/vanparys/Documents/PhD/Codes/meq.gamma';
addpath(meqpath)
[Lm,Xm] = fbt('ana',13,0,'iterq',50,'nr',66, 'nz', 64);
Ym = fbtt(Lm,Xm);

% Run EQUIL to match (might take some time)

[L,X,Y,match] = meq_to_equil(Lm,Ym,'Ns',10, 'debug',4);


c_fbt = [0.0,0.6,0.6];
c_equil = [0.9,0.4,0.3];
figure; hold on; axis equal; grid on
h_fbt = plot(NaN,NaN,'-','Color',c_fbt,'LineWidth',1.5);
h_equil = plot(NaN,NaN,'--','Color',c_equil,'LineWidth',1.5);
levels = linspace(0,match.psiN,11);
contour(Lm.rrx,Lm.zzx,(Ym.Fx-Ym.FA)/(Ym.FB-Ym.FA),levels, 'Color', c_fbt,'LineWidth',1.5)
contour(Lm.rrx,Lm.zzx,(Ym.Fx-Ym.FA)/(Ym.FB-Ym.FA),[1 1],'-.', 'Color', c_fbt,'LineWidth',1.5)
contour(match.R0*Y.RR,match.Z0+match.R0*Y.ZZ, ...
    match.psiN*Y.psiN.*ones(size(Y.RR)),levels,'--', 'Color', c_equil,'LineWidth',1.5)
legend([h_fbt,h_equil],{'MEQ','EQUIL'}, ...
    'Interpreter','latex','Box','off','FontSize',14)
xlabel('$R$ [m]','Interpreter','latex','FontSize',14)
ylabel('$Z$ [m]','Interpreter','latex','FontSize',14)
title('$\psi_N$','Interpreter','latex','FontSize',14)

% % Cool MEQ plot of the flux surfaces
% 
% [Lm,Xm] = fbt('ana',13,0);
% Ym = fbtt(Lm,Xm);
% F=Ym.Fx(2:end-1,2:end-1); r=Lm.ry; z=Lm.zy;
% C=contourc(r,z,F,[Ym.FB Ym.FB]); k=1; P=polyshape;
% while k<size(C,2), n=C(2,k); Q=regions(polyshape(C(1,k+1:k+n),C(2,k+1:k+n))); for j=1:numel(Q), if isinterior(Q(j),Ym.rA,Ym.zA), P=Q(j); break; end, end, if area(P)>0, break; end, k=k+n+1; end
% [RR,ZZ]=meshgrid(r,z); in=reshape(isinterior(P,RR(:),ZZ(:)),size(F)); climv=[min(F(in)) max(F(in))];
% 
% figure; ax=gca; hold(ax,'on'); ax.SortMethod='childorder';
% surf(Lm.rry,Lm.zzy,zeros(size(F)),F,'EdgeColor','w','LineWidth',4,'FaceColor','flat'); clim(climv)
% B=polyshape([r(1) r(end) r(end) r(1)],[z(1) z(1) z(end) z(end)]); T=triangulation(subtract(B,P));
% patch('Faces',T.ConnectivityList,'Vertices',T.Points,'FaceColor','w','EdgeColor','none')
% surf(Lm.rry,Lm.zzy,zeros(size(F)),F,'FaceColor','none','EdgeColor','k','LineWidth',.4)
% contour(Lm.rry,Lm.zzy,F,[Ym.FB Ym.FB],'k:','LineWidth',1.5)
% view(2); axis equal; box on