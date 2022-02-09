%%%% Launch inference using EP %%%%%

close all 
clear all
path(path, 'src/');
path(path, 'data/');

%% Inference using Nanchen 2006 data 
load('nanchen_2006.mat');
data = nanchen_2006.data;
Ndata = length(data);


%%% EP parameters
damp = 0.7;
max_iter = 5e4;
minvar = 1e-50;
maxvar = 1e50;
type = 'unique';
gamma_init = 1e5;
precision = 1e-4;
precision_c = 1e-5;
precision_g = 2.5e-5; % Nanchen 

min_glc = 1e10;
max_glc = -1e10;
tic
for i = 1:Ndata
     [mut(:,i), st(:,i), av(:,i), va(:,i), a(:,i), d(:,i), ...
         coeff(:,i), G, nu_e, Cov, SEP] = MetaMEP(data(i).S, ...
    data(i).b, data(i).lb, data(i).ub, data(i).mean_exp, data(i).var_exp, data(i).idx_exp, type, gamma_init,...
    precision, precision_c, precision_g, damp, max_iter, minvar, maxvar);
    %%% used of plots
    min_glc = min(min_glc, data(i).glc);
    max_glc = max(max_glc, data(i).glc);
    
end
toc

%% Plot marginal probabilities

% flux index (not measured)
idx = 59;
flux_name = data(i).rxns{idx};

%%% The colormap will mirror the value of the glucose uptake
cmap = autumn(Ndata);

plot_flux(idx, flux_name, cmap, data, mut, st, av, va, min_glc, max_glc);
 
%%
clear all

%% Inference using Schuetz 2012 data 
load('schuetz_2012.mat');
data = schuetz_2012.data;
Ndata = length(data);


%%% EP parameters
damp = 0.7;
max_iter = 5e4;
minvar = 1e-50;
maxvar = 1e50;
type = 'unique';
gamma_init = 1e5;
precision = 1e-4;
precision_c = 1e-5;
precision_g = 9e-6; % Schuetz
min_glc = 1e10;
max_glc = -1e10;
tic 
for i = 1:Ndata
     [mut(:,i), st(:,i), av(:,i), va(:,i), a(:,i), d(:,i), ...
         coeff(:,i), G, nu_e, Cov, SEP] = MetaMEP(data(i).S, ...
    data(i).b, data(i).lb, data(i).ub, data(i).mean_exp, data(i).var_exp, data(i).idx_exp, type, gamma_init, ...
    precision, precision_c, precision_g, damp, max_iter, minvar, maxvar);
    %%% used of plots
    min_glc = min(min_glc, data(i).glc);
    max_glc = max(max_glc, data(i).glc);
end
toc

%% Plot marginal probabilities

idx = 59;
flux_name = data(i).rxns{idx};
cmap = winter(Ndata);

plot_flux(idx, flux_name, cmap, data, mut, st, av, va, min_glc, max_glc);
