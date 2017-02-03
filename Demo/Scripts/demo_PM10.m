%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DYSC - DYnamic Spatiotemporal Clustering                                  %
%%%                                                                           %
%%% Authors: Francesco Finazzi and Lucia Paci                                 %
%%% E-mail: francesco.finazzi@unibg.it - lucia.paci@unicatt.it                %
%%% Affiliation: University of Bergamo - Università Cattolica del Sacro Cuore%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

addpath('../../Src');

s=RandStream('mt19937ar','Seed',5000);
RandStream.setGlobalStream(s);

load('../Data/PM10_demo.mat');

flag_plot=1; %0: do not plot outputs, 1: plot outputs

%Number of clusters
clusters=3;

%PM10 concentration
y=data(:,1:end-4);

%Station coordinates
S=data(:,end-2:-1:end-3); 

N=size(y,1); %number of stations
T=size(y,2); %number of time steps

%covariates in a NxTxB matrix
X_constant=ones(size(y,1),size(y,2)); %intercept
X_altitude=repmat(data(:,end-1),[1,T]); %station altitude
X=cat(3,X_constant,X_altitude);
B=size(X,3); %number of covariates

%dysc_data object
obj_dysc_data=dysc_data(y,X);

%dysc_site object
obj_dysc_sites=dysc_sites(S,'deg','sparse');

%dysc_hyper object
obj_dysc_hyper=dysc_hyper();

%dysc_par object
obj_dysc_par_initial=dysc_par(obj_dysc_data,clusters);
%the value of theta_zeta is 80% the maximum distance between any two sites
obj_dysc_par_initial.theta_zeta=0.8*max(obj_dysc_sites.distance_matrix(:));

%dysc_model object
obj_dysc_model=dysc_model(obj_dysc_data,obj_dysc_sites,obj_dysc_par_initial,obj_dysc_hyper);

%dysc_fit_options object
iterations=1500;
burn_in=800;
thin=1;
obj_dysc_fit_options=dysc_fit_options(iterations,burn_in,thin);

%relabelling on z (by default the relabelling is enabled)
% obj_dysc_fit_options.flag_relabelling=0;

%model parameters to be estimated (by default, all parameters are
%estimated)
% obj_dysc_fit_options.flag_update_lambda=1;
% obj_dysc_fit_options.flag_update_rho=1;
% obj_dysc_fit_options.flag_update_G=1;
% obj_dysc_fit_options.flag_update_beta=1;
% obj_dysc_fit_options.flag_update_tau=1;
% obj_dysc_fit_options.flag_update_phi=1;

%proposal type on phi
obj_dysc_fit_options.proposal_type='conditional'; %the alternative is 'random walk'
    
%tuning parameters on phi and beta when the proposal is 'rw'
if strcmp(obj_dysc_fit_options.proposal_type,'rw')
    obj_dysc_fit_options.tuning_phi=[0.005,0.005];
end
if B>0
    %the proposal on beta is random_walk so the tuning parameters are
    %needed
    obj_dysc_fit_options.tuning_beta=[0.001,0.003];
end

%model estimation
obj_dysc_model.fit(obj_dysc_fit_options);

if flag_plot
    %plot chains (burn_in and thin are optional)
    burn_in=0;
    thin=1;
    obj_dysc_model.dysc_fit_result.plot_chains(burn_in,thin);
end

%if plots are OK, aggregated results are obtained from chains
%(burn_in and thin are optional)
burn_in=0;
thin=1;
obj_dysc_model.dysc_fit_result.aggregate(burn_in,thin);

if flag_plot
    % plot densities of sigma, lambda, tau, beta, rho, G and 
    % plot the clustering result
    obj_dysc_model.dysc_fit_result.plot_density;
    obj_dysc_model.dysc_fit_result.plot_clusters;
end

%show stats
obj_dysc_model.dysc_fit_result.stats;

%compute and display DIC3 statistic
DIC3=obj_dysc_model.get_DIC3;

%save the model object for future use
save('obj_dysc_model','obj_dysc_model');

%in case you want to have more MCMC iterations you can start from the last
%iteration by accessing the dysc_par_lastiter object and by using it
%in the costructor of a new obj_dysc_model
obj_dysc_par_initial=obj_dysc_model.dysc_par_lastiter;
