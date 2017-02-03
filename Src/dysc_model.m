%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DYSC - DYnamic Spatiotemporal Clustering                                  %
%%%                                                                           %
%%% Authors: Francesco Finazzi and Lucia Paci                                 %
%%% E-mail: francesco.finazzi@unibg.it - lucia.paci@unicatt.it                %
%%% Affiliation: University of Bergamo - Università  Cattolica del Sacro Cuore%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This file is part of DYSC.
%
% DYSC is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 2 of the License, or
% (at your option) any later version.
%
% DYSC is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with DYSC. If not, see <http://www.gnu.org/licenses/>.

classdef dysc_model< handle
    
    properties (SetAccess = private)
        dysc_data=[];           %[dysc_data object]       (1x1) an object of class dysc_data
        dysc_sites=[];          %[dysc_sites object]      (1x1) an object of class dysc_sites
        dysc_par_initial=[];    %[dysc_par object]        (1x1) an object of class dysc_par which includes the initial values of all the parameters
        dysc_hyper=[];          %[dysc_hyper object]      (1x1) an object of class dysc_hyper
        
        dysc_par_simulation=[]; %[dysc_par object]        (1x1) an object of class dysc_par which includes the simulated parameters
        dysc_par_lastiter=[];   %[dysc_par object]        (1x1) an object of class dysc_par which includes the model parameters at the last MCMC iteration
        dysc_fit_result=[];     %[dysc_fit_result object] (1x1) an object of class dysc_result
        
        DIC3=[];                %[double]                 (1x1) the value of the DIC3 statistic
    
        model_estimated=0;      %[boolean]                (1x1) 0: the model is not estimated; 1: the model is estimated;
    end
    
    methods
        
        function obj = dysc_model(dysc_data,dysc_sites,dysc_par_initial,dysc_hyper)
            %DESCRIPTION: is the constructor of the class dysc_model
            %
            %INPUT
            %dysc_data:         [dysc_data object]  (1x1) an object of class dysc_data
            %dysc_sites:        [dysc_sites object] (1x1) an object of class dysc_sites
            %dysc_par_initial:  [dysc_par object]   (1x1) an object of class dysc_par which includes the initial values of all the parameters
            %<dysc_hyper>:      [dysc_hyper object] (1x1) an object of class dysc_hyper
            %
            %OUTPUT
            %obj:            [dysc_model object] (1x1) an object of class dysc_model
            
            if nargin<3
                error('The input argument dysc_data, dysc_sites, dysc_par_initial must be provided');
            end
            
            if nargin<4
                dysc_hyper=dysc_hyper();
            end
            
            if not(isa(dysc_data,'dysc_data'))
                error('dysc_data must be of class dysc_data');
            end
            if not(isa(dysc_sites,'dysc_sites'))
                error('dysc_sites must be of class dysc_sites');
            end
            if not(isa(dysc_par_initial,'dysc_par'))
                error('dysc_par must be of class dysc_par');
            end
            if not(isa(dysc_hyper,'dysc_hyper'))
                error('dysc_hyper must be of class dysc_hyper');
            end
            
            obj.dysc_data=dysc_data;
            obj.dysc_sites=dysc_sites;
            obj.dysc_par_initial=dysc_par_initial;
            obj.dysc_hyper=dysc_hyper;
        end
        
        function N = get_N(obj)
            N=obj.dysc_data.N;
        end
        
        function T = get_T(obj)
            T=obj.dysc_data.T;
        end
        
        function B = get_B(obj)
            B=obj.dysc_data.B;
        end
        
        function K = get_K(obj)
            K=obj.dysc_par_initial.K;
        end
        
        function fit(obj,dysc_fit_options)
            if not(isa(dysc_fit_options,'dysc_fit_options'))
                error('dysc_fit_options must be of class dysc_fit_options');
            end
            
            %data
            y=obj.dysc_data.y;
            X=obj.dysc_data.X;
            
            %options
            burn_in=dysc_fit_options.burn_in;
            thin=dysc_fit_options.thin;
            n_iter=dysc_fit_options.iterations;
            
            %chan_length computing
            temp=1:n_iter;
            temp=temp>burn_in&mod(temp,thin)==0;
            chain_length=sum(temp);
            
            %constants from objects
            N=obj.get_N;
            T=obj.get_T;
            B=obj.get_B;
            K=obj.get_K;
            
            %hyperparameters
            a_sigma=obj.dysc_hyper.a_sigma;
            b_sigma=obj.dysc_hyper.b_sigma;
            a_lambda=obj.dysc_hyper.a_lambda;
            b_lambda=obj.dysc_hyper.b_lambda;
            a_tau=obj.dysc_hyper.a_tau;
            b_tau=obj.dysc_hyper.b_tau;
            
            if B>0
                mu_beta0=obj.dysc_hyper.mu_beta0*ones(B,K);
                Sigma_beta0=obj.dysc_hyper.sigma_beta0*eye(B);
            end
           
            if isscalar(obj.dysc_hyper.mu_z0)
                mu_z0=obj.dysc_hyper.mu_z0*ones(K,1);
            else
                mu_z0=obj.dysc_hyper.mu_z0;
            end
            if not(size(mu_z0,1)==K&&size(mu_z0,2)==1)
                error('mu_beta0 must be Kx1');
            end
            
            if isscalar(obj.dysc_hyper.mu_z0)
                Sigma_z0=obj.dysc_hyper.sigma_z0*eye(K);
            else
                Sigma_z0=obj.dysc_hyper.sigma_z0;
            end
            if not(size(Sigma_z0,1)==K&&size(Sigma_z0,2)==K)
                error('Sigma_z0 must be KxK');
            end
            
            tuning_phi=dysc_fit_options.tuning_phi;
            if length(tuning_phi)==1
                tuning_phi=repmat(tuning_phi,[K-1,1]);
            end
            if not(length(tuning_phi)==K-1)
                error('The tuning_phi vector must be K-1x1 or a scalar');
            end
            
            tuning_beta=dysc_fit_options.tuning_beta;
            if length(tuning_beta)==1
                tuning_beta=repmat(tuning_beta,[B,1]);
            end
            if not(length(tuning_beta)==B)
                error('The tuning_beta vector must be Bx1 or a scalar');
            end
            
            %store matrices
            sigma_store=zeros(1,chain_length);
            G_store=zeros(K,chain_length);
            tau_store=zeros(K,chain_length);
            z_store=zeros(K,T+1,chain_length);
            rho_store=zeros(K-1,chain_length);
            lambda_store=zeros(K-1,chain_length);
            phi_store=zeros(N*K,T+1,chain_length);
            w_store=ones(N,T+1,chain_length);
            if B>0
                beta_store=zeros(B,K,chain_length);
            end
            theta_zeta_store=zeros(1,chain_length);
            
            %initialization of model parameters
            sigma=obj.dysc_par_initial.sigma;
            G=obj.dysc_par_initial.G;
            tau=obj.dysc_par_initial.tau;
            z=obj.dysc_par_initial.z;
            rho=obj.dysc_par_initial.rho;
            lambda=obj.dysc_par_initial.lambda;
            phi=obj.dysc_par_initial.phi;
            if strcmp(dysc_fit_options.proposal_type,'rw')
                phi(N+1:N*K,1:T+1)=unifrnd(0,1,N*(K-1),T+1);
            end
            w=obj.dysc_par_initial.w;
            beta=obj.dysc_par_initial.beta;
            theta_zeta=obj.dysc_par_initial.theta_zeta;
            if isempty(theta_zeta)
                error('The parameter theta_zeta is not updated so you must provide a value for theta_zeta in the dysc_par_initial object');
            end
            
            %acceptance rates
            acceptance_rate_phi=zeros(T+1,K-1);
            if B>0
                acceptance_rate_beta=zeros(K-1,1);
            end
            
            %constraints
            phi(1:N,:,:)=0;
            phi(:,1,:)=0;
            if B>0
                beta(:,1,:)=0;
                beta_last=beta;  
            end
            
            if (dysc_fit_options.flag_update_theta_zeta==0)
                Gamma=exp(-obj.dysc_sites.distance_matrix/theta_zeta(1,1));
                Gamma_inv=Gamma\eye(N);
            else
                error('Estimation of theta_zeta and Delta is not yet supported');
            end
            
            %hyperparameter for z1
            Sigma_z0_inv=Sigma_z0\eye(K);
            Sigma_z0_inv_mu_z0=Sigma_z0_inv*mu_z0;
            
            %hyperparameter for phi1
            phi0=obj.dysc_hyper.mu_phi0*ones(N,1);
            
            %Cholesky of Gamma is computed only one time. This
            %works only if theta_zeta is fixed and not estimated.
            chol_Gamma=cholcov(Gamma,0);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %           MCMC            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('MCMC started...')
            aTn2_sigma =a_sigma+T*N/2;
            aTn2_lambda=a_lambda+T*N/2;
            aT2_tau    =a_tau+(T-1)/2;
            
            %matrix for z selection
            wsel=ones(N,1)*(0:T-1)*K;
            n_z=zeros(K,1);
            y_z=zeros(K,1);
            
            acceptance_phi_vector=zeros(500,K-1);
            acceptance_beta_vector=zeros(500,K-1);
            acceptance_phi_index=1;
            acceptance_beta_index=1;
            last_line_length=0;
            prog_bar_length=50;
            step_iter=3;
            step_iter_vec=[1 2 5 10 20 25 50 100 150 200 250 500 1000 2000];
            refresh_period=2; %in seconds
            time1=clock;
            
            %%%%%%%%%%%% mcmc started %%%%%%%%%%%%%
            counter_store=1;
            for i=1:n_iter
                %feedback to user
                if mod(i,step_iter)==0||i==n_iter
                    time2=clock;
                    elapsed_time=etime(time2,time1);
                    time_per_iter=elapsed_time/i;
                    iter_per_sec=(i-2)/elapsed_time;
                    step_iter=ceil(iter_per_sec)*refresh_period;
                    [~,idx]=min(abs(step_iter_vec-step_iter));
                    step_iter=step_iter_vec(idx);
                    time_to_finish=round((n_iter-i)*time_per_iter);
                    completed=i/n_iter;
                    
                    if time_to_finish<60
                        time_to_finish_string=[num2str(time_to_finish,2),'s'];
                    else
                        if time_to_finish<3600
                            m=floor(time_to_finish/60);
                            s=time_to_finish-60*m;
                            time_to_finish_string=[num2str(m,2),'m ',num2str(s,2),'s'];
                        else
                            h=floor(time_to_finish/3600);
                            time_to_finish=time_to_finish-h*3600;
                            m=floor(time_to_finish/60);
                            s=time_to_finish-60*m;
                            time_to_finish_string=[num2str(h,2),'h ',num2str(m,2),'m ',num2str(s,2),'s'];
                        end
                    end
                    fprintf(1, repmat('\b',1,last_line_length));
                    to_display=['Iteration ',num2str(i),'/',num2str(n_iter),' (',num2str(time_per_iter,2),'s per iter - ',time_to_finish_string,' to finish)'];
                    segments=round(completed*prog_bar_length);
                    perc=[num2str(round(completed*100)),'%'];
                    perc_length=length(perc);
                    if segments+perc_length<=prog_bar_length
                        to_display=[to_display,sprintf('\n'),'[',repmat('>',1,segments),perc,repmat('-',1,prog_bar_length-segments-perc_length),']'];
                    else
                        to_display=[to_display,sprintf('\n'),'[',repmat('>',1,segments),repmat('-',1,prog_bar_length-segments),']'];
                    end
                    if dysc_fit_options.flag_update_phi==1
                        for k=2:K
                            instant_rate=sum(acceptance_phi_vector(1:acceptance_phi_index-1,k-1))/acceptance_phi_index;
                            if instant_rate==0
                                instant_rate_string='0.0';
                            else
                                if instant_rate==1
                                    instant_rate_string='1.0';
                                else
                                    instant_rate_string=num2str(instant_rate,2);
                                    if strcmp(instant_rate_string,'0')
                                        instant_rate_string='0.0';
                                    else if strcmp(instant_rate_string,'1')
                                            instant_rate_string='1.0';
                                        else
                                            while length(instant_rate_string)<4
                                                instant_rate_string=[instant_rate_string,'0'];
                                            end
                                        end
                                    end
                                end
                            end
                            gauge_pos=round(instant_rate*prog_bar_length);
                            to_display=[to_display,sprintf('\n'),'Acceptance rate on phi_',num2str(k),sprintf('\n'),'[',repmat('-',1,gauge_pos),char(9829),repmat('-',1,prog_bar_length-gauge_pos-1),']'];
                            to_display=[to_display,sprintf('\n'),repmat(' ',1,gauge_pos-1),instant_rate_string];
                        end
                    end
                    if dysc_fit_options.flag_update_beta==1
                        to_display=[to_display,sprintf('\n')];
                        for k=2:K
                            instant_rate=sum(acceptance_beta_vector(1:acceptance_beta_index-1,k-1))/acceptance_beta_index;
                            if instant_rate==0
                                instant_rate_string='0.0';
                            else
                                if instant_rate==1
                                    instant_rate_string='1.0';
                                else
                                    instant_rate_string=num2str(instant_rate,2);
                                    if strcmp(instant_rate_string,'0')
                                        instant_rate_string='0.0';
                                    else if strcmp(instant_rate_string,'1')
                                            instant_rate_string='1.0';
                                        else
                                            while length(instant_rate_string)<4
                                                instant_rate_string=[instant_rate_string,'0'];
                                            end
                                        end
                                    end
                                end
                            end
                            gauge_pos=round(instant_rate*prog_bar_length);
                            to_display=[to_display,sprintf('\n'),'Acceptance rate on beta_{k=',num2str(k),'}',sprintf('\n'),'[',repmat('-',1,gauge_pos),char(9829),repmat('-',1,prog_bar_length-gauge_pos-1),']'];
                            to_display=[to_display,sprintf('\n'),repmat(' ',1,gauge_pos-1),instant_rate_string];
                        end
                    end
                    
                    disp(to_display);
                    last_line_length=length(to_display)+1;
                end
                
                %sigma
                if dysc_fit_options.flag_update_sigma
                    ww=w(:,2:end)+wsel;
                    zz=z(:,2:end);
                    res_sigma=y-zz(ww);
                    sigma=1/gamrnd(aTn2_sigma,1/(b_sigma+0.5*sum(diag(res_sigma'*res_sigma))));
                end
                
                %lambda
                if dysc_fit_options.flag_update_lambda
                    for k=2:K
                        res_lambda=phi((k-1)*N+1:k*N,2:end)-rho(k-1,1)*phi((k-1)*N+1:k*N,1:end-1);
                        lambda(k-1,1)=1/gamrnd(aTn2_lambda,1/(b_lambda+0.5*sum(diag(res_lambda'*Gamma_inv*res_lambda))));
                        %note that Gamma_inv is used in the above line as Gamma_inv is common to the K-1 \phi
                    end
                end
                
                %tau
                if dysc_fit_options.flag_update_tau
                    res_tau=z(:,3:end)-repmat(G(:,1),[1,T-1]).*z(:,2:end-1); %modificato
                    tau(:,1)=1./gamrnd(aT2_tau,1./(b_tau+0.5*sum(res_tau.^2,2)));
                end
                
                %z
                if dysc_fit_options.flag_update_z
                    GSigmaG=diag((1./tau(:,1)).*diag((G(:,1)*G(:,1)')));
                    for t=2:T+1
                        w_ti=w(:,t);
                        for k=1:K
                            n_z(k)=sum(w_ti==k);
                            L=w_ti==k;
                            y_z(k)=sum(y(L,t-1));
                        end
                        
                        if t==T+1
                            V_z=diag(diag(1./(diag(n_z)/sigma+diag(1./tau(:,1)))));
                            M_z=y_z/sigma+G(:,1).*(1./tau(:,1)).*z(:,t-1);
                        else
                            if t==2
                                V_z=diag(diag(1./(diag(n_z)/sigma+GSigmaG+Sigma_z0_inv)));
                                M_z=y_z/sigma+G(:,1).*(1./tau(:,1)).*z(:,t+1)+Sigma_z0_inv_mu_z0;
                            else
                                V_z=diag(diag(1./(diag(n_z)/sigma+GSigmaG+diag(1./tau(:,1)))));
                                M_z=y_z/sigma+G(:,1).*(1./tau(:,1)).*z(:,t-1)+G(:,1).*(1./tau(:,1)).*z(:,t+1);
                            end
                        end
                        
                        z(:,t)=mvnrnd(V_z*M_z,V_z);
                        if dysc_fit_options.flag_relabelling==1
                            %ordering constraint
                            z(:,t)=sort(z(:,t));
                        end
                    end
                end
                
                %G
                if dysc_fit_options.flag_update_G
                    v_G=1./(sum(z(:,3:end).^2,2)./tau(:,1)+1/obj.dysc_hyper.sigma_G);
                    m_G=sum(z(:,2:end-1).*z(:,3:end),2)./tau(:,1);
                    
                    for k=1:K
                        pd = makedist('Normal',v_G(k)*m_G(k),sqrt(v_G(k)));
                        try
                            pdt = truncate(pd,-1,1);
                            G(k,1)=random(pdt,1);
                        catch
                            if v_G(k)*m_G(k)>1
                                G(k,1)=1-eps;
                            else
                                G(k,1)=-1+eps;
                            end
                        end
                    end
                end
                
                %rho
                if dysc_fit_options.flag_update_rho
                    V_rho=zeros(K-1);
                    M_rho=zeros(K-1,1);
                    for t=2:T+1
                        for k=2:K
                            chi=phi((k-1)*N+1:k*N,t-1);
                            temp=chi'*Gamma_inv/lambda(k-1,1);
                            V_rho(k-1,k-1)=V_rho(k-1,k-1)+temp*chi;
                            M_rho(k-1)=M_rho(k-1)+temp*phi((k-1)*N+1:k*N,t);
                        end
                    end
                    V_rho=(V_rho+1/obj.dysc_hyper.sigma_rho.*eye(K-1))\eye(K-1);
                    rho(:,1)=mvnrnd((V_rho*M_rho),V_rho,1);
                    
                    if not(isequal(diag(diag(V_rho)),V_rho))
                        error('V_rho must be diagonal to use univariate truncated normals');
                    end
                    for k=1:K-1
                        pd = makedist('Normal',V_rho(k,k)*M_rho(k),sqrt(V_rho(k,k)));
                        try
                            pdt = truncate(pd,-1,1);
                            rho(k,1)=random(pdt,1);
                        catch
                            if V_rho(k,k)*M_rho(k)>1
                                rho(k,1)=1-eps;
                            else
                                rho(k,1)=-1+eps;
                            end
                        end
                    end
                end
                
                %phi metropolis
                if dysc_fit_options.flag_update_phi
                    if strcmp(dysc_fit_options.proposal_type,'rw')
                        for t=2:T+1
                            H=zeros(N,K);
                            for h=1:K
                                L=w(:,t)==h;
                                H(L,h)=1;
                            end
                            
                            %This is needed as input of mnpdf_fast
                            L_H=H==1;
                            
                            Phi=reshape(phi(:,t),[N,K]);%PHI(:,:,t);
                            Phi_past=reshape(phi(:,t-1),[N,K]);%PHI(:,:,t-1);
                            if t<T+1
                                Phi_fut=reshape(phi(:,t+1),[N,K]);%PHI(:,:,t+1);
                            end
                            
                            if acceptance_phi_index>500
                                acceptance_phi_vector=circshift(acceptance_phi_vector,[-1 0]);
                                acceptance_phi_index=500;
                            end
                            
                            % updating for each k=2,...,K
                            for k=2:K
                                %chol_lambda_Gamma is the Cholesky of lambda*Gamma
                                chol_lambda_Gamma=sqrt(lambda(k-1,1))*chol_Gamma;
                                logSqrtDet_chol = sum(log(diag(chol_lambda_Gamma)));
                                
                                % current value
                                if B>0
                                    X_Phi=squeeze(X(:,t-1,:))*beta+Phi;
                                else
                                    X_Phi=Phi;
                                end
                                temp=exp(X_Phi);
                                Pi=temp./repmat(sum(temp,2),[1,K]);
                                
                                % candidate
                                % compute quantities of the conditional prior
                                chol_tuning_Gamma=sqrt(tuning_phi(k-1))*chol_Gamma;
                                Phi_new=Phi;
                                Phi_new(:,k)=(randn(1,size(chol_lambda_Gamma,1))* chol_tuning_Gamma)'+Phi(:,k);
                                
                                if B>0
                                    X_Phi_new=squeeze(X(:,t-1,:))*beta+Phi_new;
                                else
                                    X_Phi_new=Phi_new;
                                end
                                temp=exp(X_Phi_new);
                                Pi_new=temp./repmat(sum(temp,2),[1,K]);
                                
                                f_phi=rho(k-1,1)*Phi_past(:,k);
                                %log rate
                                if t==T+1
                                    l_r=log(mvnpdf_chol(Phi_new(:,k),f_phi,chol_lambda_Gamma,logSqrtDet_chol))-... %p(\phi_new|\phi_t-1)
                                        log(mvnpdf_chol(Phi(:,k),f_phi,chol_lambda_Gamma,logSqrtDet_chol))+...$p(\phi_corrente|\phi_t-1)
                                        sum(log(mnpdf_fast(L_H,Pi_new)))-...
                                        sum(log(mnpdf_fast(L_H,Pi)));
                                else
                                    if t==2
                                        l_r=log(mvnpdf_chol(Phi_new(:,k),phi0,chol_lambda_Gamma,logSqrtDet_chol))-... %p(\phi_new|\phi_t-1)
                                            log(mvnpdf_chol(Phi(:,k),phi0,chol_lambda_Gamma,logSqrtDet_chol))+...$p(\phi_corrente|\phi_t-1)
                                            log(mvnpdf_chol(Phi_fut(:,k),rho(k-1,1)*Phi_new(:,k),chol_lambda_Gamma,logSqrtDet_chol))-... %p(\phi_new|\phi_t-1)
                                            log(mvnpdf_chol(Phi_fut(:,k),rho(k-1,1)*Phi(:,k),chol_lambda_Gamma,logSqrtDet_chol))+...$p(\phi_corrente|\phi_t-1)
                                            sum(log(mnpdf_fast(L_H,Pi_new)))-...
                                            sum(log(mnpdf_fast(L_H,Pi)));
                                    else
                                        l_r=log(mvnpdf_chol(Phi_new(:,k),f_phi,chol_lambda_Gamma,logSqrtDet_chol))-... %p(\phi_new|\phi_t-1)
                                            log(mvnpdf_chol(Phi(:,k),f_phi,chol_lambda_Gamma,logSqrtDet_chol))+...$p(\phi_corrente|\phi_t-1)
                                            log(mvnpdf_chol(Phi_fut(:,k),rho(k-1,1)*Phi_new(:,k),chol_lambda_Gamma,logSqrtDet_chol))-... %p(\phi_new|\phi_t-1)
                                            log(mvnpdf_chol(Phi_fut(:,k),rho(k-1,1)*Phi(:,k),chol_lambda_Gamma,logSqrtDet_chol))+...$p(\phi_corrente|\phi_t-1)
                                            sum(log(mnpdf_fast(L_H,Pi_new)))-...
                                            sum(log(mnpdf_fast(L_H,Pi)));
                                    end
                                end
                                
                                if log(rand)<l_r
                                    %accept
                                    Phi(:,k)=Phi_new(:,k);
                                    
                                    acceptance_rate_phi(t,k-1)=acceptance_rate_phi(t,k-1)+1;
                                    acceptance_phi_vector(acceptance_phi_index,k-1)=1;
                                else
                                    acceptance_phi_vector(acceptance_phi_index,k-1)=0;
                                end
                            end
                            acceptance_phi_index=acceptance_phi_index+1;
                            phi(:,t)=reshape(Phi,[N*K,1]);
                        end
                    else
                        %VVV_chol_temp is computed only once as invariant with respect to t
                        VVV_chol_temp=cell(K-1,1);
                        for k=2:K
                            VVV_chol_temp{k-1}=sqrt(lambda(k-1,1)/(1+rho(k-1,1)^2))*chol_Gamma;
                        end
                        
                        for t=2:T+1
                            H=zeros(N,K);
                            for h=1:K
                                L=w(:,t)==h;
                                H(L,h)=1;
                            end
                            
                            %This is needed as input of mnpdf_fast
                            L_H=H==1;
                            
                            Phi=reshape(phi(:,t),[N,K]);%PHI(:,:,t);
                            Phi_past=reshape(phi(:,t-1),[N,K]);%PHI(:,:,t-1);
                            if t<T+1
                                Phi_fut=reshape(phi(:,t+1),[N,K]);%PHI(:,:,t+1);
                            end
                            
                            % updating for each k=2,...,K
                            acceptance_in_cycle=1;
                            old_numerator_computed=0;
                            old_denominator_computed=0;
                            
                            if acceptance_phi_index>500
                                acceptance_phi_vector=circshift(acceptance_phi_vector,[-1 0]);
                                acceptance_phi_index=500;
                            end
                            
                            for k=2:K
                                %This block must stay in the cycle 2:K asi Phi is possibly
                                %updated in the cycle. However, computation can be avoided
                                %if Phi is not updated in the cycle.
                                
                                %current value
                                if acceptance_in_cycle==1
                                    if B>0
                                        X_Phi=squeeze(X(:,t-1,:))*beta+Phi;
                                    else
                                        X_Phi=Phi;
                                    end
                                    temp=exp(X_Phi);
                                    Pi=temp./repmat(sum(temp,2),[1,K]);
                                end
                                
                                % candidate
                                if t==T+1
                                    %When t=T+1, VVV_chol is VVV_chol_temp not divided by (1+rho(k,i)^2);
                                    VVV_chol=VVV_chol_temp{k-1}*sqrt(1+rho(k-1,1)^2);
                                    MMM=rho(k-1,1)*Phi_past(:,k);
                                else
                                    VVV_chol=VVV_chol_temp{k-1};
                                    if t==2
                                        MMM=rho(k-1,1)/(rho(k-1,1)^2+1)*(Phi_fut(:,k));
                                    else
                                        MMM=rho(k-1,1)/(rho(k-1,1)^2+1)*(Phi_fut(:,k)+Phi_past(:,k));
                                    end
                                end
                                
                                Phi_new=Phi;
                                Phi_new(:,k)=(randn(1,size(VVV_chol,1))*VVV_chol)'+MMM;
                                
                                if B>0
                                    X_Phi_new=squeeze(X(:,t-1,:))*beta+Phi_new;
                                else
                                    X_Phi_new=Phi_new;
                                end
                                temp=exp(X_Phi_new);
                                Pi_new=temp./repmat(sum(temp,2),[1,K]);
                                
                                %log rate
                                num=sum(log(mnpdf_fast(L_H,Pi_new)));
                                if acceptance_in_cycle==1&&old_numerator_computed==1
                                    den=old_num;
                                else
                                    if acceptance_in_cycle==0&&old_denominator_computed==1
                                        den=old_den;
                                    else
                                        den=sum(log(mnpdf_fast(L_H,Pi)));
                                    end
                                end
                                l_r=num-den;
                                
                                if log(rand)<l_r
                                    %accept
                                    Phi(:,k)=Phi_new(:,k);
                                    acceptance_phi_vector(acceptance_phi_index,k-1)=1;
                                    acceptance_rate_phi(t,k-1)=acceptance_rate_phi(t,k-1)+1;
                                    %the next denominator is the current nominator
                                    old_num=num;
                                    acceptance_in_cycle=1;
                                    old_numerator_computed=1;
                                else
                                    %the next denominator is the current denominator
                                    old_den=den;
                                    acceptance_phi_vector(acceptance_phi_index,k-1)=0;
                                    acceptance_in_cycle=0;
                                    old_denominator_computed=1;
                                end
                            end
                            acceptance_phi_index=acceptance_phi_index+1;
                            phi(:,t)=reshape(Phi,[N*K,1]);
                        end
                    end
                end
                
                %w
                if dysc_fit_options.flag_update_w
                    for t=2:T+1
                        Phi=reshape(phi(:,t),[N,K]);
                        if B>0
                            X_Phi=squeeze(X(:,t-1,:))*beta+Phi;
                        else
                            X_Phi=Phi;
                        end
                        temp=exp(X_Phi);
                        Pi=temp./repmat(sum(temp,2),[1,K]);
                        
                        %likelihood to be in cluster k
                        like_final=normpdf(repmat(y(:,t-1),[1,K]),repmat(z(:,t)',[N,1]),repmat(sqrt(sigma),[N,K]));
                        
                        Pi_star_num=Pi.*like_final;
                        Pi_star=Pi_star_num./repmat(sum(Pi_star_num,2),[1,K]);
                        
                        temp=cumsum(Pi_star,2);
                        r=rand(N,1);
                        H=zeros(N,K);
                        for k=1:K
                            H(:,k)=temp(:,k)>=r;
                        end
                        
                        HH=H';
                        [~,idx] = max(HH(:,sum(HH)>0));
                        w(:,t)=idx';
                    end
                end
                
                %beta
                if B>0
                    if dysc_fit_options.flag_update_beta
                        XX=reshape(X,[N*T,B]);
                        Phi_b=zeros(N*T,K);
                        for k=1:K
                            temp=phi((k-1)*N+1:k*N,2:end);
                            Phi_b(:,k)=temp(:);
                        end
                        
                        HH=zeros(N,T,K);
                        for t=2:T+1
                            for k=1:K
                                L=w(:,t)==k;
                                HH(L,t-1,k)=1;
                            end
                        end
                        HH=reshape(HH,[N*T,K]);
                        L_H=HH==1;
                        
                        beta_temp=beta_last;
                        for k=2:K
                            beta_temp(:,k)=mvnrnd(beta_last(:,k),diag(tuning_beta))';
                            
                            num=exp(XX*beta_last+Phi_b);
                            den=repmat(sum(exp(XX*beta_last+Phi_b),2),[1,K]);
                            Pi=num./den;
                            
                            num_new=exp(XX*beta_temp+Phi_b);
                            den_new=repmat(sum(exp(XX*beta_temp+Phi_b),2),[1,K]);
                            Pi_new=num_new./den_new;
                            
                            %Can be optimized as it was done for phi
                            log_r=sum(log(mnpdf_fast(L_H,Pi_new)))-...
                                sum(log(mnpdf_fast(L_H,Pi)))+...
                                log(mvnpdf(beta_temp(:,k),mu_beta0(:,k),Sigma_beta0))-...
                                log(mvnpdf(beta_last(:,k),mu_beta0(:,k),Sigma_beta0));
                            
                            if acceptance_beta_index>500
                                acceptance_beta_vector(:,k-1)=circshift(acceptance_beta_vector(:,k-1),[-1 0]);
                                acceptance_beta_index=500;
                            end
                            
                            if log(rand)>log_r
                                %reject
                                beta_temp(:,k)=beta_last(:,k);
                                acceptance_beta_vector(acceptance_beta_index,k-1)=0;
                            else
                                acceptance_rate_beta(k-1)=acceptance_rate_beta(k-1)+1;
                                acceptance_beta_vector(acceptance_beta_index,k-1)=1;
                            end
                        end
                        acceptance_beta_index=acceptance_beta_index+1;
                        beta=beta_temp;
                        beta_last=beta;
                    end
                end
                
                %theta_zeta
                if dysc_fit_options.flag_update_theta_zeta
                    error('At the moment the software does not allow to estimate theta_zeta');
                end
                
                if mod(i,thin)==0&&i>burn_in
                    sigma_store(:,counter_store)=sigma;
                    lambda_store(:,counter_store)=lambda;
                    tau_store(:,counter_store)=tau;
                    z_store(:,:,counter_store)=z;
                    G_store(:,counter_store)=G;
                    rho_store(:,counter_store)=rho;
                    phi_store(:,:,counter_store)=phi;
                    w_store(:,:,counter_store)=w;
                    if B>0
                        beta_store(:,:,counter_store)=beta;
                    end
                    theta_zeta_store(:,counter_store)=theta_zeta;
                    counter_store=counter_store+1;
                end
            end
            disp('MCMC ended.')
            
            obj.dysc_fit_result=dysc_fit_result(obj.dysc_sites);
            obj.dysc_fit_result.sigma=sigma_store;
            obj.dysc_fit_result.G=G_store;
            obj.dysc_fit_result.tau=tau_store;
            obj.dysc_fit_result.rho=rho_store;
            obj.dysc_fit_result.lambda=lambda_store;
            obj.dysc_fit_result.z=z_store;
            obj.dysc_fit_result.phi=phi_store;
            obj.dysc_fit_result.w=w_store;
            if B>0
                obj.dysc_fit_result.beta=beta_store;
            end
            obj.dysc_fit_result.theta_zeta=theta_zeta_store;
            
            %create an object for the parameters at the last MCMC iteration
            obj.dysc_par_lastiter=dysc_par(obj.dysc_data,K);
            obj.dysc_par_lastiter.sigma=sigma_store(:,end);
            obj.dysc_par_lastiter.G=G_store(:,end);
            obj.dysc_par_lastiter.tau=tau_store(:,end);
            obj.dysc_par_lastiter.rho=rho_store(:,end);
            obj.dysc_par_lastiter.lambda=lambda_store(:,end);
            obj.dysc_par_lastiter.z=z_store(:,:,end);
            obj.dysc_par_lastiter.phi=phi_store(:,:,end);
            obj.dysc_par_lastiter.w=w_store(:,:,end);
            if B>0
                obj.dysc_par_lastiter.beta=beta_store(:,:,end);
            end
            obj.dysc_par_lastiter.theta_zeta=theta_zeta_store(:,end);
            
            %compute acceptance rate
            obj.dysc_fit_result.acceptance_rate_phi=acceptance_rate_phi./n_iter;
            if B>0
                obj.dysc_fit_result.acceptance_rate_beta=acceptance_rate_beta/n_iter;
            end
            obj.model_estimated=1;
        end
        
        function DIC3 = get_DIC3(obj)
            %compute the DIC3 statistic
            
            if obj.model_estimated==0
                error('The model is not estimated. Call the method fit of the dysc_model object');
            end
            if obj.dysc_fit_result.result_aggregated==0
                error('Model fittin results are not aggregated. Call the aggregate method of the dysc_fit_result object first');
            end

            N=obj.get_N;
            T=obj.get_T;
            B=obj.get_B;
            K=obj.get_K;
            nsample=obj.dysc_fit_result.get_chain_length;
            
            log_mix=zeros(N,T,nsample);
            f_hat=zeros(N,T);
            temp=zeros(K,nsample);
            
            for t=2:T+1
                Phi=reshape(obj.dysc_fit_result.phi(:,t,:),[N,K,nsample]);
                for s=1:N
                    if B>0
                        for k=2:K
                            temp(k,:)=squeeze(obj.dysc_data.X(s,t-1,:))'*squeeze(obj.dysc_fit_result.beta(:,k,:));
                        end
                        X_Phi=temp+squeeze(Phi(s,:,:));
                    else
                        X_Phi=squeeze(Phi(s,:,:));
                    end
                    num=exp(X_Phi);
                    Pi=num./repmat(sum(num,1),[K,1]);
                    
                    %Normal density
                    density=normpdf(repmat(obj.dysc_data.y(s,t-1),[K,nsample]),squeeze(obj.dysc_fit_result.z(:,t,:)),repmat(sqrt(obj.dysc_fit_result.sigma),[K,1]));
                    
                    % sum pi x N
                    SpN=sum(Pi.*density,1);
                    
                    % log_mixture
                    log_mix(s,t-1,:)=log(SpN);
                    
                    %f_hat
                    f_hat(s,t-1)=mean(SpN);
                end
            end
            
            Dbar=-2/nsample*sum(sum(sum(log_mix)));
            obj.DIC3=2*Dbar+2*(sum(sum(log(f_hat))));
            DIC3=obj.DIC3;
            disp(['DIC3: ',num2str(DIC3)]);
        end
        
        function sim(obj)
            %Simulate data using the parameters in the dysc_par_initial
            %object. The result is stored in the y property of the
            %dysc_data object.
            par=obj.dysc_par_initial;
            
            N=obj.get_N;
            T=obj.get_T;
            K=par.K;
            
            %phi
            phi0=obj.dysc_hyper.mu_phi0*ones(N,1);
            Phi=zeros(N,K,T+1);
            for k=2:K
                Gamma=par.lambda(k-1)*exp(-obj.dysc_sites.distance_matrix/par.theta_zeta);
                Phi(:,k,1)=phi0;
                for t=2:T+1
                    Phi(:,k,t)=par.rho(k-1)*Phi(:,k,t-1)+mvnrnd(zeros(N,1),Gamma)';
                end
            end
            
            %pi
            Pi=zeros(N,K,T);
            for t=2:T+1
                Pi(:,:,t-1)=exp(Phi(:,:,t))./repmat(sum(exp(Phi(:,:,t)),2),[1,K]);
            end
            
            %H
            H=zeros(N,K,T);
            for t=1:T
                H(:,:,t)=mnrnd(ones(N,1),Pi(:,:,t));
            end
            
            %Z
            z=zeros(K,T+1);
            z(:,1)=obj.dysc_hyper.mu_z0;
            for t=2:T+1
               z(:,t)=par.G.*z(:,t-1)+mvnrnd(zeros(K,1),diag(par.tau))';
            end
            %A constant is added to separate the components of Z
            constant=(1:K)*2;
            for t=1:T+1
                z(:,t)=z(:,t)+constant';
            end
            
            %Y
            y=zeros(N,T);
            w=zeros(N,T+1);
            for t=1:T
                HH=H(:,:,t)';
                [~,wt] = max(HH(:,sum(HH)>0));
                y(:,t)=z(wt,t+1)+normrnd(0,sqrt(par.sigma),N,1);
                w(:,t+1)=wt;
            end
            obj.dysc_data.y=y; 
            
            obj.dysc_par_simulation=obj.dysc_par_initial;
            obj.dysc_par_simulation.phi=reshape(Phi,N*K,T+1);
            obj.dysc_par_simulation.z=z;
            obj.dysc_par_simulation.w=w;
        end   
    end
end

        