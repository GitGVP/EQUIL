function basis = assemble_variational_bspline_axis_regular(m,nq,p,Ns,Nh)
% B-splines with the regular polar Taylor class enforced on the axis span.

if m < 2 || m ~= round(m) || p < max([3,Ns,Nh])
    error('Use m >= 2 and spline_p >= max(3,Ns,Nh).');
end

knots = [zeros(1,p+1),(1:m-1)/m,ones(1,p+1)];
nb = m+p;
[g,w] = lgwt(nq,0,1);
g = flip(g(:));
w = flip(w(:));
r = zeros(m*nq,1);
wr = zeros(m*nq,1);
for element = 1:m
    rows = (element-1)*nq+(1:nq);
    left = (element-1)/m;
    r(rows) = left+g/m;
    wr(rows) = w/m;
end

[N0,N1] = bspline_eval_all(knots,p,r);
[E0,E1] = bspline_eval_all(knots,p,[0,1]);
axis_derivatives = first_span_derivatives(knots,p,m,nb);

nprofiles = 3+Ns+Nh;
% The linear P coefficient is the radial-normalization gauge.  Its
% geometrical remainder has the regular first-harmonic powers r^3,r^5,...
leading_power = [2,2,1,1:Ns,1:Nh];
basis.B0 = cell(nprofiles,1);
basis.B1 = cell(nprofiles,1);
basis.B0_axis = cell(nprofiles,1);
basis.B1_axis = cell(nprofiles,1);
basis.B0_edge = cell(nprofiles,1);
basis.B1_edge = cell(nprofiles,1);
basis.lift0 = cell(Ns+Nh,1);
basis.lift1 = cell(Ns+Nh,1);
basis.lift0_axis = zeros(Ns+Nh,1);
basis.lift1_axis = zeros(Ns+Nh,1);
basis.lift0_edge = ones(Ns+Nh,1);
basis.lift1_edge = zeros(Ns+Nh,1);
basis.profile_lengths = zeros(1,nprofiles);

edge_lift = zeros(nb,1);
edge_lift(end) = 1;
for profile = 1:nprofiles
    fixed_edge = profile > 3;
    transform = constrained_transform( ...
        axis_derivatives,p,nb,leading_power(profile),fixed_edge);
    basis.B0{profile} = sparse(N0.'*transform);
    basis.B1{profile} = sparse(N1.'*transform);
    basis.B0_axis{profile} = sparse(E0(:,1).'*transform);
    basis.B1_axis{profile} = sparse(E1(:,1).'*transform);
    basis.B0_edge{profile} = sparse(E0(:,2).'*transform);
    basis.B1_edge{profile} = sparse(E1(:,2).'*transform);
    basis.profile_lengths(profile) = size(transform,2);
    if fixed_edge
        shape = profile-3;
        basis.lift0{shape} = sparse(N0.'*edge_lift);
        basis.lift1{shape} = sparse(N1.'*edge_lift);
        basis.lift0_axis(shape) = E0(:,1).'*edge_lift;
        basis.lift1_axis(shape) = E1(:,1).'*edge_lift;
        basis.lift0_edge(shape) = E0(:,2).'*edge_lift;
        basis.lift1_edge(shape) = E1(:,2).'*edge_lift;
    end
end

basis.r = r;
basis.w = wr;
basis.profile_starts = cumsum([1,basis.profile_lengths(1:end-1)]);
end

function transform = constrained_transform(D,p,nb,n,fixed_edge)
% Keep only r^n,r^(n+2),... Taylor powers on the axis span.
orders = [0:n-1,n+1:2:p];
local = null(D(orders+1,1:p+1));
last_free = nb-fixed_edge;
tail = (p+2):last_free;
transform = zeros(nb,size(local,2)+numel(tail));
transform(1:p+1,1:size(local,2)) = local;
if ~isempty(tail)
    transform(tail,size(local,2)+(1:numel(tail))) = eye(numel(tail));
end
end

function D = first_span_derivatives(knots,p,m,nb)
% Polynomial coefficients in u=m*r give exact zero-derivative constraints.
u = (0:p)'/p;
r = u/m;
N = bspline_eval_all(knots,p,r);
V = ones(p+1,p+1);
for power = 1:p
    V(:,power+1) = u.^power;
end
coefficients = V\N.';
D = zeros(p+1,nb);
for order = 0:p
    D(order+1,:) = factorial(order)*coefficients(order+1,:);
end
end
