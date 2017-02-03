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

classdef dysc_fit_result < handle
    
    %constant
    %I: number of MCMC iterations
    
    properties
        %chains
        sigma     =[];   %[double]      (QxI)
        G         =[];   %[double]      (KxI)
        tau       =[];   %[double]      (KxI)
        rho       =[];   %[double]      (K-1xI)
        lambda    =[];   %[double]      (K-1xI)
        z         =[];   %[double]      (KxT+1xI)
        phi       =[];   %[double]      (N*KxT+1xI)
        w         =[];   %[double]      (NxT+1xI)
        beta      =[];   %[double]      (BxKxI)
        theta_zeta=[];   %[double>0]    (1xI)
        
        acceptance_rate_phi=[];     %[double] (Tx1)
        acceptance_rate_beta=[];    %[double] (Kx1)
    end
    
    properties (SetAccess = private)
        %averages and modes
        sigma_avg     =[];   %[double]      (Qx1)
        G_avg         =[];   %[double]      (Kx1)
        tau_avg       =[];   %[double]      (Kx1)
        rho_avg       =[];   %[double]      (K-1x1)
        lambda_avg    =[];   %[double]      (K-1x1)
        z_025         =[];   %[double]      (KxT+1)
        z_500         =[];   %[double]      (KxT+1)
        z_975         =[];   %[double]      (KxT+1)
        z_avg         =[];   %[double]      (KxT+1)
        phi_avg       =[];   %[double]      (N*KxT+1)
        w_mode        =[];   %[double]      (NxT+1)
        beta_avg      =[];   %[double]      (BxK)
        theta_zeta_avg=[];   %[double>0]    (1x1)
        
        w_var         =[];   %[double]      (NxT+1) variability of the w
        
        result_aggregated=0; %[integer]     (1x1) 0: aggregation not done, 1: aggregation done
    end
    
    properties (SetAccess = private, GetAccess = private)
        dysc_sites=[];       %[dysc_sites object]  (1x1) a dysc_sites object
    end
    
    methods
        
        function obj = dysc_fit_result(dysc_sites)
            %DESCRIPTION: is the constructor of the class dysc_fit_result
            %
            %INPUT
            %dysc_sites:        [dysc_sites object]       (1x1) an object of class dysc_sites
            %
            %OUTPUT
            %obj:               [dysc_fit_result object]  (1x1) an object of class dysc_hyper
            if not(isa(dysc_sites,'dysc_sites'))
                error('dysc_sites must be of class dysc_sites');
            end
            obj.dysc_sites=dysc_sites;
        end
        
        function n = get_chain_length(obj)
            %get the length of the chains
            n=length(obj.sigma);
        end
        
        function aggregate(obj,burn_in,thin)
            %DESCRIPTION: compute summaries from the chains
            %
            %INPUT
            %burn_in:   [integer>0] (1x1) the number of burn-in iterations
            %thin:      [integer>0] (1x1) thinning
            %OUTPUT
            %obj:       [dysc_fit_result object]  (1x1) an object of class dysc_hyper
            if nargin<2
                burn_in=0;
            end
            if nargin<3
                thin=1;
            end
            
            disp('Aggregation started...');
            obj.sigma=obj.sigma(:,burn_in+1:thin:end);
            obj.G=obj.G(:,burn_in+1:thin:end);
            obj.tau=obj.tau(:,burn_in+1:thin:end);
            obj.rho=obj.rho(:,burn_in+1:thin:end);
            obj.lambda=obj.lambda(:,burn_in+1:thin:end);
            obj.z=obj.z(:,:,burn_in+1:thin:end);
            obj.phi=obj.phi(:,:,burn_in+1:thin:end);
            obj.w=obj.w(:,:,burn_in+1:thin:end);
            if not(isempty(obj.beta))
                obj.beta=obj.beta(:,:,burn_in+1:thin:end);
            end
            obj.theta_zeta=obj.theta_zeta(:,burn_in+1:thin:end);
            
            %aggregation
            obj.sigma_avg = mean(obj.sigma,2);
            obj.G_avg = mean(obj.G,2);
            obj.tau_avg = mean(obj.tau,2);
            obj.rho_avg = mean(obj.rho,2);
            obj.lambda_avg = mean(obj.lambda,2);
            %z
            obj.z_avg = mean(obj.z,3);
            obj.z_025=quantile(obj.z,0.025,3);
            obj.z_500=quantile(obj.z,0.500,3);
            obj.z_975=quantile(obj.z,0.975,3);
            %phi
            obj.phi_avg = mean(obj.phi,3);
            %w
            obj.w_mode = mode(obj.w,3);
            obj.w_var=zeros(size(obj.w,1),size(obj.w,2));
            for i=1:size(obj.w,1)
                for j=1:size(obj.w,2)
                    tbl=tabulate(squeeze(obj.w(i,j,:)));
                    if size(tbl,1)>1
                        mode_freq=max(tbl(:,2));
                        max_val=sum(tbl(:,2))*(size(tbl,1)-1);
                        obj.w_var(i,j)=sum(mode_freq-tbl(:,2))/max_val;
                    else
                        obj.w_var(i,j)=1;
                    end
                end
            end
            %beta
            if not(isempty(obj.beta))
                obj.beta_avg = mean(obj.beta,3);
            end
            %theta_zeta
            obj.theta_zeta_avg = mean(obj.theta_zeta,2);
            
            obj.result_aggregated=1;
            disp('Aggregation ended.');
        end
        
        function plot_chains(obj,burn_in,thin,lag)
            %DESCRIPTION: plot the chains
            %
            %INPUT
            %obj:       [dysc_fit_result]   (1x1)
            %burn_in:   [integer>0]         (1x1) the number of burn-in iterations
            %thin:      [integer>0]         (1x1) thinning
            %lag:       [integer>0]         (1x1) number of autocorrelation lags in graphs
            %OUTPUT
            %chain plots
            
            if nargin<2
                burn_in=0;
            end
            if nargin<3
                thin=1;
            end
            if nargin<4
                lag=20;
            end
            
            K=size(obj.lambda,1)+1;
            B=size(obj.beta,1);
            T=size(obj.z,2)-1;
            I=size(obj.lambda,2);
            
            L=burn_in+1:thin:I;
            
            %beta
            if B>0
                for b=1:B
                    figure('Name',['beta_',num2str(b),' chains'],'NumberTitle','off');
                    for k=1:K-1
                        subplot(2,K-1,k);
                        plot(squeeze(obj.beta(b,k+1,L)));
                        title(['beta_{b=',num2str(b),',k=',num2str(k+1),'} - burn in: ',num2str(burn_in),', thin: ',num2str(thin)]);
                        xlabel('Sample');
                        grid on
                        xlim([1,length(L)]);
                        subplot(2,K-1,k+K-1);
                        autocorr(squeeze(obj.beta(b,k+1,L)),lag);
                        title(['Autocor. of beta_{b=',num2str(b),',k=',num2str(k+1),'} - burn in: ',num2str(burn_in),', thin: ',num2str(thin)]);
                        xlabel('Lag');
                        grid on
                    end
                end
            end
            
            % plot rho
            figure('Name','rho chains','NumberTitle','off');
            for k=1:K-1
                subplot(2,K-1,k);
                plot(obj.rho(k,L));
                title(['rho_',num2str(k+1),' - burn in: ',num2str(burn_in),', thin: ',num2str(thin)]);
                xlabel('Sample');
                grid on
                xlim([1,length(L)]);
                subplot(2,K-1,k+K-1);
                autocorr(obj.rho(k,L),lag);
                title(['Autocor. of rho_',num2str(k+1),' - burn in: ',num2str(burn_in),', thin: ',num2str(thin)]);
                xlabel('Lag');
                grid on
            end
            
            %lambda
            figure('Name','lambda chains','NumberTitle','off');
            for k=1:K-1
                subplot(2,K-1,k);
                plot(obj.lambda(k,L));
                title(['lambda_',num2str(k+1),' - burn in: ',num2str(burn_in),', thin: ',num2str(thin)]);
                xlabel('Sample');
                grid on
                xlim([1,length(L)]);
                subplot(2,K-1,k+K-1);
                autocorr(obj.lambda(k,L),lag);
                title(['Autocor. of lambda_',num2str(k+1),' - burn in: ',num2str(burn_in),', thin: ',num2str(thin)]);
                xlabel('Lag');
                grid on
            end
           
            %G
            figure('Name','G chains','NumberTitle','off');
            for k=1:K
                subplot(2,K,k);
                plot(obj.G(k,L));
                title(['g_',num2str(k),' - burn in: ',num2str(burn_in),', thin: ',num2str(thin)]);
                xlabel('Sample');
                grid on
                xlim([1,length(L)]);
                subplot(2,K,k+K);
                autocorr(obj.G(k,L),lag);
                title(['Autocor. of g_',num2str(k),' - burn in: ',num2str(burn_in),', thin: ',num2str(thin)]);
                xlabel('Lag');
                grid on
            end
            
            %tau
            figure('Name','tau chains','NumberTitle','off');
            for k=1:K
                subplot(2,K,k);
                plot(obj.tau(k,L));
                title(['tau_',num2str(k),' - burn in: ',num2str(burn_in),', thin: ',num2str(thin)]);
                xlabel('Iteration');
                grid on
                xlim([1,length(L)]);
                subplot(2,K,k+K);
                autocorr(obj.tau(k,L),lag);
                title(['Autocor. of tau_',num2str(k),' - burn in: ',num2str(burn_in),', thin: ',num2str(thin)]);
                xlabel('Lag');
                grid on
            end
            
            %sigma
            figure('Name','sigma chain','NumberTitle','off');
            subplot(1,2,1);
            plot(obj.sigma(L));
            title(['sigma - burn in: ',num2str(burn_in),', thin: ',num2str(thin)]);
            xlabel('Sample');
            grid on
            xlim([1,length(L)]);
            subplot(1,2,2);
            autocorr(obj.sigma(burn_in+1:thin:end),lag);
            title(['Autocor. of sigma - burn in: ',num2str(burn_in),', thin: ',num2str(thin)]);
            xlabel('Lag');
            grid on
            
            % plot z
            w_g=4;
            h_g=3;
            pages=ceil(T/(w_g*h_g));
            for p=1:pages
                figure('Name','z chains','NumberTitle','off');
                for t=(p-1)*w_g*h_g+2:min(p*w_g*h_g,T)+1
                    subplot(h_g,w_g,t-1-(p-1)*w_g*h_g);
                    plot(squeeze(obj.z(:,t,L))');
                    title(['z_{t=',num2str(t-1),'} - burn in: ',num2str(burn_in),', thin: ',num2str(thin)]);
                    xlim([1,length(L)]);
                end
            end
            
            disp('Plots of phi and w are not automatically provided. You can plot them by accessing the phi and w properties of this object');
        end
        
        function plot_clusters(obj)
            %DESCRIPTION: plot the clustering result
            %
            %INPUT
            %obj:       [dysc_fit_result]   (1x1)
            %OUTPUT
            %clustering result
            
            if obj.result_aggregated==1
                %plot of z
                figure('Name','z clustering','NumberTitle','off');
                hold on
                for i=1:size(obj.z_avg,1) 
                    switch i
                        case 1
                            style='d';
                            color='g';
                        case 2
                            style='d';
                            color='b';
                        case 3
                            style='d';
                            color='r';
                        case 4
                            style='d';
                            color='c';
                        case 5
                            style='d';
                            color='m';
                        case 6
                            style='d';
                            color='y';
                    end
                    plot(obj.z_avg(i,2:end)',[style,color]); 
                    width=size(obj.z_avg,2)/300;
                    for t=2:size(obj.z_avg,2)
                        plot([t-1 t-1],[obj.z_025(i,t),obj.z_975(i,t)],['-',color]);
                        plot([t-1-width t-1+width],[obj.z_025(i,t),obj.z_025(i,t)],['-',color]);
                        plot([t-1-width t-1+width],[obj.z_975(i,t),obj.z_975(i,t)],['-',color]);
                    end
                end
                title('Latent temporal components z (95% CI)');
                xlabel('Time step');
                xlim([0,size(obj.z_avg,2)]);
                grid on
                
                %plot of w
                figure('Name','w clustering','NumberTitle','off');
                T=size(obj.w_mode,2)-1;
                r=ceil(sqrt(9/16*T));
                c=ceil(T/r);
                
                z_min=min(obj.z_avg(:));
                z_max=max(obj.z_avg(:));
                for t=2:T+1
                    subplot(r,c,t-1);
                    w_temp=obj.w_mode(:,t);
                    w_alpha=obj.w_var(:,t);
                    if strcmp(obj.dysc_sites.type,'grid')
                        w_alpha=reshape(w_alpha,obj.dysc_sites.grid_size);
                        w_temp=reshape(w_temp,obj.dysc_sites.grid_size);
                        imagesc(w_temp,'AlphaData',w_alpha);
                        set(gca,'Ydir','Normal');
                        hold on
                        x_c=size(w_temp,2)+0.5;
                        y_c=size(w_temp,1);
                        for i=1:size(obj.z_avg(:,t),1)
                            switch i
                                case 1
                                    style='o';
                                    color='g';
                                case 2
                                    style='d';
                                    color='b';
                                case 3
                                    style='*';
                                    color='r';
                                case 4
                                    style='^';
                                    color='c';
                                case 5
                                    style='p';
                                    color='m';
                                case 6
                                    style='p';
                                    color='y';
                            end
                            plot(x_c,(obj.z_avg(i,t)-z_min)/(z_max-z_min)*y_c+0.5,style,'MarkerFaceColor',color,'MarkerEdgeColor',color)
                        end        
                    else
                        for i=1:max(w_temp)
                            switch i
                                case 1
                                    style='.g';
                                case 2
                                    style='.b';
                                case 3
                                    style='.r';
                                case 4
                                    style='.c';
                                case 5
                                    style='.m';
                                case 6
                                    style='.y';
                            end
                            L=w_temp==i;
                            plot(obj.dysc_sites.coordinates(L,1),obj.dysc_sites.coordinates(L,2),style);
                            hold on
                        end
                    end
                    title(['w(s,t=',num2str(t-1),')']);
                    axis equal
                    axis tight
                end  
            else
                error('No results available for plot. Call the aggregate method of the class discy_fit_result first.');
            end
        end
        
        function plot_density(obj)
            %DESCRIPTION: plot the densities
            %
            %INPUT
            %obj:       [dysc_fit_result]   (1x1)
            %OUTPUT
            %none
            
            if obj.result_aggregated==1
                K=size(obj.lambda,1)+1;
                B=size(obj.beta,1);
                
                %beta
                if B>0
                    for b=1:B
                        figure('Name',['beta_',num2str(b),' density'],'NumberTitle','off');
                        for k=1:K-1
                            subplot(1,K-1,k);
                            hist(squeeze(obj.beta(b,k+1,:)),40);
                            title(['Histogram of beta_{b=',num2str(b),',k=',num2str(k+1),'} - Mean: ',num2str(obj.beta_avg(b,k))]);
                        end
                    end
                end
                
                %rho
                figure('Name','rho density','NumberTitle','off');
                for k=1:K-1
                    subplot(1,K-1,k);
                    hist(obj.rho(k,:),40);
                    title(['Histogram of rho_',num2str(k+1),' - Mean: ',num2str(obj.rho_avg(k))]);
                end
                
                %lambda
                figure
                for k=1:K-1
                    subplot(1,K-1,k);
                    hist(obj.lambda(k,:),40);
                    title(['Histogram of lambda_',num2str(k+1),' - Mean: ',num2str(obj.lambda_avg(k))]);
                end
                
                
                %G
                figure('Name','G density','NumberTitle','off');
                for k=1:K
                    subplot(1,K,k);
                    hist(obj.G(k,:),40);
                    title(['Histogram of g_',num2str(k),' - Mean: ',num2str(obj.G_avg(k))]);
                end
                
                %tau
                figure('Name','tau density','NumberTitle','off');
                for k=1:K
                    subplot(1,K,k);
                    hist(obj.tau(k,:),40);
                    title(['Histogram of tau_',num2str(k),' - Mean: ',num2str(obj.tau_avg(k))]);
                end
                
                %sigma
                figure('Name','sigma density','NumberTitle','off');
                hist(obj.sigma,40);
                title(['Histogram of sigma. Mean: ',num2str(obj.sigma_avg)]);
            else
                error('No results available for plot. Call the aggregate method of the class dysc_fit_result first.');
            end
        end
        
        function tab = stats(obj)
            %DESCRIPTION: provide stats for the estimated parameters
            %
            %INPUT
            %obj:       [dysc_fit_result]   (1x1)
            %OUTPUT
            %tab:       [matlab table]      (1x1)
            
            if obj.result_aggregated==1
                sample_size=length(obj.lambda);
                K=size(obj.lambda,1)+1;
                B=size(obj.beta,1);
                
                counter_tab=1;
                %beta
                if B>0
                    for b=1:B
                        for k=1:K-1
                            tab(counter_tab,1)=mean(squeeze(obj.beta(b,k+1,:)));
                            tab(counter_tab,2)=std(squeeze(obj.beta(b,k+1,:)));
                            tab(counter_tab,3:5)=quantile(squeeze(obj.beta(b,k+1,:)),[0.025 0.500 0.975]);
                            tab(counter_tab,6)=sample_size;
                            label{counter_tab,1}=['beta_{b=',num2str(b),',k=',num2str(k+1),'}'];
                            counter_tab=counter_tab+1;
                        end
                    end
                end
                %rho
                for k=1:K-1
                    tab(counter_tab,1)=mean(obj.rho(k,:));
                    tab(counter_tab,2)=std(obj.rho(k,:));
                    tab(counter_tab,3:5)=quantile(obj.rho(k,:),[0.025 0.500 0.975]);
                    tab(counter_tab,6)=sample_size;
                    label{counter_tab,1}=['rho_',num2str(k+1)];
                    counter_tab=counter_tab+1;
                end
                %lambda
                for k=1:K-1
                    tab(counter_tab,1)=mean(obj.lambda(k,:));
                    tab(counter_tab,2)=std(obj.lambda(k,:));
                    tab(counter_tab,3:5)=quantile(obj.lambda(k,:),[0.025 0.500 0.975]);
                    tab(counter_tab,6)=sample_size;
                    label{counter_tab,1}=['lambda_',num2str(k+1)];
                    counter_tab=counter_tab+1;
                end
                %G
                for k=1:K
                    tab(counter_tab,1)=mean(obj.G(k,:));
                    tab(counter_tab,2)=std(obj.G(k,:));
                    tab(counter_tab,3:5)=quantile(obj.G(k,:),[0.025 0.500 0.975]);
                    tab(counter_tab,6)=sample_size;
                    label{counter_tab,1}=['g_',num2str(k)];
                    counter_tab=counter_tab+1;
                end
                %tau
                for k=1:K
                    tab(counter_tab,1)=mean(obj.tau(k,:));
                    tab(counter_tab,2)=std(obj.tau(k,:));
                    tab(counter_tab,3:5)=quantile(obj.tau(k,:),[0.025 0.500 0.975]);
                    tab(counter_tab,6)=sample_size;
                    label{counter_tab,1}=['tau_',num2str(k)];
                    counter_tab=counter_tab+1;
                end
                %sigma
                tab(counter_tab,1)=mean(obj.sigma);
                tab(counter_tab,2)=std(obj.sigma);
                tab(counter_tab,3:5)=quantile(obj.sigma,[0.025 0.500 0.975]);
                tab(counter_tab,6)=sample_size;
                label{counter_tab,1}='sigma';
                counter_tab=counter_tab+1;
                %theta_zeta
                tab(counter_tab,1)=mean(obj.theta_zeta);
                tab(counter_tab,2)=std(obj.theta_zeta);
                tab(counter_tab,3:5)=quantile(obj.theta_zeta,[0.025 0.500 0.975]);
                tab(counter_tab,6)=sample_size;
                label{counter_tab,1}='theta_zeta';
                
                tab=round(tab*1000)/1000;
                tab=array2table(tab);
                tab.Properties.VariableNames={'Mean','Std','p025','p500','p975','sample'};
                tab.Properties.RowNames=label;
                tab
            else
                error('Stats are not available. Call the aggregate method of the class dysc_fit_result first.');
            end
        end
    end
end