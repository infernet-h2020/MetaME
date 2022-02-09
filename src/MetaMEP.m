function [mut, st, av, va, a, d, c_all, G, nu_e, Cov, SEP] = MetaMEP(S, ...
    b, nuinf, nusup, av_exp_in, va_exp_in, idx_exp, type, givengamma, ...
    precision, precision_c, precision_g, damp, max_iter, minvar, maxvar)



Nr = size(S,2);
Nm = size(S,1);
Nexp = length(idx_exp);
% Full set of EP parameters
a = zeros(Nr,1,'double');
d = ones(Nr,1,'double');
% \nu moments
av = zeros(Nr,1,'double');
va = ones(Nr,1,'double');
mut = zeros(Nr,1,'double');
st = ones(Nr,1,'double');
% scaling
factor = max(max(abs(nuinf)), max(abs(nusup)));
nusup = nusup / factor;
nuinf = nuinf / factor;
b = b / factor;
damp_c = damp;

% row echelon form
[aux, idxd] = rref(cat(2, S, b));
Nd = length(idxd);
aux = aux(1:Nd,:);
idxf = setdiff(1:Nr,idxd);
Nf = length(idxf);
A = aux(:,idxf);
Ap = A';
y = aux(:,end);
basis = zeros(Nr,Nf);
basis(idxd,:) = -aux(:,idxf);
basis(idxf,:) = eye(Nf, Nf);


% parameters of the approximation for free and dependent fluxes
Df = sparse(1:Nf, 1:Nf, 1./d(idxf));
Dd = sparse(1:Nd, 1:Nd, 1./d(idxd));

% indices of free and dependent fluxes among all fluxes
posd = ones(Nr,1) * NaN;
posf = ones(Nr,1) * NaN;
for i = 1:Nr
    if ismember(i,idxd)
        posd(i) = find(idxd == i);
    else
        posf(i) = find(idxf == i);
    end
end

% Lagrange multipliers - means (with holes)
cf = zeros(Nf, 1,'double');
cd = zeros(Nd, 1,'double');


% Split measured fluxes in free and dep.
Nf_e = nnz(ismember(idxf, idx_exp));
Nd_e = nnz(ismember(idxd, idx_exp));
Ad_exp = zeros(Nd_e, Nf);
yd_exp = zeros(Nd_e, 1);
av_exp_d = zeros(Nd_e, 1);
va_exp_d = zeros(Nd_e, 1);
va_exp_f = zeros(Nf_e, 1);
av_exp_f = zeros(Nf_e, 1);
% indices of measured fluxes within the vectors idxd and idxf
idxf_e = zeros(Nf_e,1); 
idxd_e = zeros(Nd_e,1);

k = 1;
for i = 1:Nd
    if ismember(idxd(i), idx_exp)
       aux = find(idx_exp == idxd(i));
       Ad_exp(k,:) = A(i,:);
       yd_exp(k) = y(i);
       av_exp_d(k) = av_exp_in(aux) / factor;
       va_exp_d(k) = va_exp_in(aux) / (factor^2);
       idxd_e(k) = i;
       k = k + 1;
    end
end
Adp_exp = Ad_exp';
k = 1;
for i = 1:Nf
    if ismember(idxf(i),idx_exp)
        aux = find(idx_exp == idxf(i));
        av_exp_f(k) = av_exp_in(aux) / factor;
        va_exp_f(k) = va_exp_in(aux) / (factor^2);
        idxf_e(k) = i;
        k = k + 1;
    end
end
gamma0 = givengamma;
% Lagrange multipliers - variances 
Gd = givengamma .* eye(Nd_e);
Gf = givengamma .* eye(Nf_e);
Gf_full = eye(Nf);
diag_f = zeros(Nf,1);
diag_f(idxf_e) = givengamma;
Gf_full = Gf_full - diag(diag(Gf_full)) + diag(diag_f);

dampg = 0.1;
iter = 0;
err = 100;
err_g = 100;
err_c = 100;
err_lin = 100;
while ((err > precision || err_c > precision_c || err_g > precision_g) && iter < max_iter)
    
    iter = iter + 1;
    
    % update EP parameters, fixed cf and cd
    Sigma1 = Df + Ap * Dd * A;
    Sigma = inv(Sigma1);
    s = diag(Sigma);
    mu = Sigma1 \ (Df * a(idxf) + Ap * Dd * (y - a(idxd)) + cf - Ap*cd);
    oldvar = va;
    oldav = av;
    
    for i = 1:Nr
        if(ismember(i,idxf))
            idx = posf(i);
            ss = min(maxvar, max(minvar, Sigma(idx,idx)));
            muv = mu(idx);
        else	
            idx = posd(i);
            x = -A(idx,:)';
            ss = min(maxvar, max(minvar, x'*Sigma*x));
            muv = x'*mu + y(idx);
        end
        st(i) = min(maxvar, max(minvar, 1/(1/ss - 1/d(i))));
        mut(i) = (muv-(a(i)*ss)/d(i))/(1-ss/d(i));
        if(ss == d(i))
    		mut(i) = 0.5*(nuinf(i) + nusup(i));
        end
        
        s05(i) = st(i)^0.5;
        x0 = (nuinf(i)-mut(i))/s05(i); 
        x1 = (nusup(i)-mut(i))/s05(i);     
        [z,eps] = compute_mom_vec(x0, x1);    
        av(i) = mut(i) + z * s05(i);
        va(i) = max(0, st(i) * (1 + eps));
    end
    
    err = max(max(abs(av-oldav)), max(abs(va-oldvar)));
    err_lin = 1/Nr * (S*av - b)'*(S*av - b);
   
    
    % update cd, cf
    vL = inv(Df + Ap * Dd * A + Adp_exp * Gd * Ad_exp + Gf_full);
    
    U = Gd - Gd * Ad_exp * vL * Adp_exp * Gd;
    V = Gd * Ad_exp * vL(:,idxf_e) * Gf;
    Y = Gf - Gf * vL(idxf_e,idxf_e) * Gf;
    
    M1 = [U V; V' Y];
    u_d = Gd * yd_exp - Gd * Ad_exp * vL * (Df * a(idxf) + Ap * Dd * (y - a(idxd)) + Adp_exp * Gd * yd_exp);
    u_f = Gf * vL(idxf_e, :) * (Df *a(idxf) + Ap * Dd * (y - a(idxd)) + Adp_exp * Gd * yd_exp);
    
       
    res = M1*cat(1,av_exp_d, av_exp_f) - cat(1,u_d, u_f); % new coeff
    
    m_all = M1 \ ( cat(1,cd(idxd_e),cf(idxf_e))  +  cat(1,u_d, u_f) ); % first moments, model
    md = m_all(1:Nd_e);
    mf = m_all(Nd_e+1:end);
    
    err_c = mean(abs(m_all - cat(1, av_exp_d, av_exp_f)));
    % coeff update
    cd(idxd_e) = damp_c * cd(idxd_e) + (1-damp_c) * res(1:Nd_e);
    cf(idxf_e) = damp_c * cf(idxf_e) + (1-damp_c) * res(Nd_e+1:end);
    
   % pause
    % moments matching (update "a" and "d")
    new_d = 1./(1./va - 1./st);
    new_d = min(maxvar, max(minvar, new_d));
    new_a = av + (av - mut) .* new_d.* (1./st);    

    a = damp * a + (1-damp) * new_a;
    d = damp * d + (1-damp) * new_d;
    Df = sparse(1:Nf, 1:Nf, 1./d(idxf));
    Dd = sparse(1:Nd, 1:Nd, 1./d(idxd));
    
    if err_c < precision_c
        % update gamma
        if strcmp(type, 'fixed')
            err_g = 0;
        else
            T1 = eye(Nf_e) - vL(idxf_e, idxf_e) * Gf ...
                - vL(idxf_e,:) * Adp_exp * Gd * inv(Gd - Gd*Ad_exp*vL*Adp_exp*Gd) * Gd * Ad_exp * vL(:,idxf_e) * Gf;
            T = inv(T1);
            W1 = eye(Nd_e) - Ad_exp * vL * Adp_exp * Gd ...
                - Ad_exp * vL(:,idxf_e) * Gf * inv(Gf - Gf * vL(idxf_e,idxf_e) *Gf) * Gf * vL(idxf_e,:) * Adp_exp * Gd;
            W = inv(W1);
            gd = diag(W) ./ va_exp_d;
            gf = diag(T) ./ va_exp_f;

            gamma0 = mean([sum(diag(W)) / sum(va_exp_d), sum(diag(T)) / sum(va_exp_f)]);
            va_e_all = diag(inv(M1));
            if strcmp(type, 'unique') 
                    Gd = Gd - diag(diag(Gd)) + dampg * diag(diag(Gd)) + (1-dampg) * gamma0 * eye(Nd_e);
                    Gf = Gf - diag(diag(Gf)) + dampg * diag(diag(Gf)) + (1-dampg) * gamma0 * eye(Nf_e);
                    diag_f(idxf_e) = gamma0;
                    Gf_full = Gf_full - diag(diag(Gf_full)) + dampg * diag(diag(Gf_full)) + (1-dampg) * diag(diag_f);
                    err_g = abs( mean(va_e_all) - mean(cat(1, va_exp_d, va_exp_f)) );

            else
                    Gd = Gd - diag(diag(Gd)) + dampg * diag(diag(Gd)) + (1-dampg) * diag(gd);
                    Gf = Gf - diag(diag(Gf)) + dampg * diag(diag(Gf)) + (1-dampg) * diag(gf);
                    diag_f(idxf_e) = gf;
                    Gf_full = Gf_full - diag(diag(Gf_full)) + dampg * diag(diag(Gf_full)) + (1-dampg) * diag(diag_f);
                    err_g = max( abs(va_e_all - cat(1, va_exp_d, va_exp_f) ));
            end
           
            fprintf('it:%i err_lin:%1.2e conv_EP:%1.2e err_c:%1.2e err_g:%1.2e gamma:%1.2e rank:%i\n',...
                iter, err_lin, err, err_c, err_g, gamma0,rank(M1));
           
        end
    end
    if ~mod(iter,100)
        fprintf('it:%i err_lin:%1.2e conv_EP:%1.2e err_c:%1.2e err_g:%1.2e gamma:%1.2e\n', iter, err_lin, err, err_c, err_g, gamma0);
    end
end
    fprintf('it:end err_lin:%1.2e conv_EP:%1.2e err_c:%1.2e err_g:%1.2e gamma:%1.2e\n', err_lin, err, err_c, err_g, gamma0);

   
    % \nu statistics
    st = st * factor^2;
    mut = mut * factor;
    a = a*factor;
    d = d*factor^2;
    av = av*factor;
    va = va*factor^2;
    b = b * factor;
    
    Df = sparse(1:Nf, 1:Nf, 1./d(idxf));
    Dd = sparse(1:Nd, 1:Nd, 1./d(idxd));
    Sigma1 = Df + Ap * Dd * A;
    Covf = inv(Sigma1);
    Cov = basis * Covf * basis';
    orth_matrix = orth_norm(find(idxf == size(S,2)),basis);
    SEP = entropy(orth_matrix' * Cov * orth_matrix);
    
    % Lagrange multipliers
    c_all = zeros(Nr,1);
    G = sparse(Nr,Nr);
    
    for i = 1:Nr
        if (ismember(i, idxf(idxf_e)))
            aux = find(idxf == i);
            c_all(i) = cf(aux) / factor;
            G(i,i) = Gf_full(aux,aux) / (factor^2);
        end
        if (ismember(i, idxd(idxd_e)))
            aux = find(idxd == i);
            aux1 = find(idxd(idxd_e) == i);
            c_all(i) = cd(aux) / factor;
            G(i,i) = Gd(aux1,aux1) / (factor^2);
        end
    end
    
    %\nu^e statistics
    
    nu_e = zeros(Nexp,1);
    Cov_e = zeros(Nexp,1);
    cov_aux = diag(inv(M1));
    for i = 1:Nexp
        if ismember(idx_exp(i), idxd(idxd_e))
            aux = find(idxd(idxd_e) == idx_exp(i));
            nu_e(i) = m_all(aux) * factor;
            Cov_e(i) = cov_aux(aux) * factor * factor;
        end
        if ismember(idx_exp(i), idxf(idxf_e))
            aux = find(idxf(idxf_e) == idx_exp(i));
            nu_e(i) = m_all(Nd_e + aux) * factor;
            Cov_e(i) = cov_aux(Nd_e + aux) * factor * factor;
        end
    end
end

function S = entropy(M)

    N = size(M,1);
    [L,p] = chol(M);
    if(~p)
        S = sum(log(diag(L))) + 0.5*N*log(2*pi*exp(1));
    else
        S = NaN;
    end
end

function y = Phi(x)
    y = 0.5 * (1 + erf(x/sqrt(2)));
end

function y = phi(x)
    y = 1/sqrt(2*pi) * exp(-x.^2/2);
end

function [z, z1] = compute_mom_vec(xinf, xsup)
    n = size(xinf, 1);
    z = zeros(n,1,'double');
    z1 = zeros(n,1,'double');
    for i = 1:n
        [z(i),z1(i)]=compute_mom5d(xinf(i),xsup(i));
    end
end

function [scra1, scra12] = compute_mom5d(xinf, xsup)
    if xsup - xinf < 1e-10
        scra1 = 0.5*(xsup + xinf);
        scra12 = -1;
        return;
    end
    
    if min(abs(xinf), abs(xsup)) <= 6. || xinf*xsup <= 0
        Phisup   = Phi(xsup);
        phisup   = phi(xsup);
        Phiinf   = Phi(xinf);
        phiinf   = phi(xinf);
        scra1 = (phiinf - phisup)/(Phisup - Phiinf);
        scra2 = (xinf * phiinf - xsup*phisup)/(Phisup - Phiinf);        
    else
        delta2 = (xsup^2 - xinf^2)*0.5;
        if delta2 > 40.
            scra1 = xinf^5/(3 - xinf^2 + xinf^4);
            scra2 = xinf^6/(3 - xinf^2 + xinf^4);
        else
            scra1 = (xinf*xsup)^5 * (1. - exp(delta2)) / (-exp(delta2)*(3.0-xinf^2 + xinf^4)*xsup^5 + xinf^5*(3-xsup^2 + xsup^4));
            scra2 = (xinf*xsup)^5 * (xsup - xinf*exp(delta2)) / (-exp(delta2)*(3.0-xinf^2 + xinf^4)*xsup^5 + xinf^5*(3-xsup^2 + xsup^4));
        end
    end
    scra12=scra2-scra1^2;
end

       




