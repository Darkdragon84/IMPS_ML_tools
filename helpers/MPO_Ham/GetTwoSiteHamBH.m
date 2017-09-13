function [H,bop] = GetTwoSiteHamBH(d,param)

% gives the nearest neighbor Hamiltonian for the Bose Hubbard model

% hopping
if isfield(param,'t') && ~isempty(param.t) && abs(param.t)>0,t=param.t;
else t=-1;
end

% onsite pair interaction
if isfield(param,'U') && ~isempty(param.U) && abs(param.U)>0,U=param.U;
else error('U not specified');
end

% NN density density interaction
if isfield(param,'V') && ~isempty(param.V) && abs(param.V)>0,V=param.V;
else V=0;
end

% chemical potential
if isfield(param,'mu') && ~isempty(param.mu) && abs(param.mu)>0,mu=param.mu;
else mu=0;
end

% pair creation/annihilation term
if isfield(param,'nu') && ~isempty(param.nu) && abs(param.nu)>0,nu=param.nu;
else nu=0;
end


ndiag = sqrt(1:d-1);

bop = diag(ndiag,1);
nop = diag(0:d-1);
Id = eye(d);


onsite = U/2*nop*(nop-Id) - mu*nop + nu*(bop'*bop' + bop*bop);

H = -t*(kron(bop,bop') + kron(bop',bop)) + V*kron(nop,nop) + 0.5*(kron(onsite,Id) + kron(Id,onsite));

end

