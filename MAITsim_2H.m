function [isoStore_delta] = MAITsim_2H(Store_mat,Store_init,Flux_mat,atm_mat,time1,time2)
% Model Agnostic Isotope Tracer simulator v1.0
%   Simulates isotope tracers (2H here) for input storages and fluxes
%   Currently works for normal fluxes, evaporation and precipitation, in
%   flux-state models.
%   
%   The model is based on conservation of isotope mass, and solved based on
%   the backwards Euler method (for numerical stability). All fluxes
%   have the same concentration as the storage they are sourced from, with
%   two current exceptions:
%   - Precipitation has the concentration of the specified atmosphere
%   - Evaporation has a concentration estimated based on the Craig-Gordon
%     equation
%
%   Model inputs (more details in Isotope model requirements doc):
%   - Store_mat: matrix of information on water storage states for all model time steps
%   - Store_init: matrix of initial isotope and storage values
%   - Flux_mat: matrix of information on fluxes for all model time steps
%   - atm_mat: matrix of information on the atmosphere for all model time steps
%   - time1: vector of start times for each time step
%   - time2: vector of end times for each time step

% Minimum stored water mass to run isotope tracer mixing:
mass_check_value=0.000001;

% Determine total counts for the model:
storage_count=size(Store_mat,1);
flux_count=size(Flux_mat,1);
time_count=length(time1);

% Step 0: initilize the storage
Store1=Store_init(:,1);
isoStore1=Store1.*delta2conc_2H(Store_init(:,2));
isoStore_conc=zeros(storage_count,time_count+1);
isoStore_conc(:,1)=delta2conc_2H(Store_init(:,2));

% Main program: simulate for each time interval
for t=1:time_count
    % Step 1: Reset for new time interval
    % Zero internal variables
    Inflow=zeros(storage_count,1);
    isoInflow=zeros(storage_count,1);
    Outflow_norm=zeros(storage_count,1);
    Outflow_evap=zeros(storage_count,1);
    isoStore2=zeros(storage_count,1);
    K_count=zeros(storage_count,1);
    J_count=zeros(storage_count,1);
    K=zeros(storage_count,flux_count);
    J=zeros(storage_count,flux_count);
    
    % Find the length of the time interval
    dt=time2(t)-time1(t);
    
    % Step 2: Identify contributing and destination storage sets
    for n=1:flux_count
        % Count and track the inflows into each storage based on the values
        %   in the flux destinations and check if flux source is modeled
        destination_i=Flux_mat(n,3,t);
        if destination_i>0
            K_count(destination_i)=K_count(destination_i)+1;
            K(destination_i,K_count(destination_i))=n;
        end
        % Count and track the outflows from each storage based on the
        %   values in the flux sources and check if flux source is modeled
        source_i=Flux_mat(n,2,t);
        if source_i>0
            J_count(source_i)=J_count(source_i)+1;
            J(source_i,J_count(source_i))=n;
        end
    end

    % Simulate each storage sequentially
    for i=1:storage_count
            % Zero internal variables
            E_ET_split=0;

        % Step 3: Add up the inputs for storage i
        for n=1:K_count(i)
            flux_id=K(i,n);
            %Total water added:
            Inflow(i)=Inflow(i)+Flux_mat(flux_id,1,t)*dt;
            if Flux_mat(flux_id,2,t)>0 
                %Modeled water source:
                %Concentration calculated at the end of this interval
                Inflow_conc=isoStore_conc(Flux_mat(flux_id,2,t),t+1);
                %Heavy isotope water added:
                isoInflow(i)=isoInflow(i)+Inflow_conc*Flux_mat(flux_id,1,t)*dt;
            else 
                %Specified source, precipitation only currently:
                %Precipitation concentration is from atm matrix
                Inflow_conc=delta2conc_2H(atm_mat(Store_mat(Flux_mat(flux_id,3,t),2,t),3,t));
                %Heavy isotope water added:
                isoInflow(i)=isoInflow(i)+Inflow_conc*Flux_mat(flux_id,1,t)*dt;
            end
        end
        
        % Step 4: Add up the outputs from storage i
        for n=1:J_count(i)
            flux_id=J(i,n);
            %Total water removed through any process except evaporation:
            Outflow_norm(i)=Outflow_norm(i)+Flux_mat(flux_id,1,t)*(1-Flux_mat(flux_id,4,t))*dt;
            %Total water removed through evaporation:
            Outflow_evap(i)=Outflow_evap(i)+Flux_mat(flux_id,1,t)*(Flux_mat(flux_id,4,t))*dt;

            % Weighted E/ET split estimate for storage i
            %   will be divided by total E if actual used for turbulence
            %   factor estimate in craig_gordon
            E_ET_split=E_ET_split+Flux_mat(flux_id,4,t)*Flux_mat(flux_id,1,t)*(Flux_mat(flux_id,4,t))*dt;
            % Warning if E/ET fraction is impossible:
            if Flux_mat(flux_id,4,t)>1 || Flux_mat(flux_id,4,t)<0
                disp(['WARNING: E/ET fraction error at time step: ',num2str(t), ', flux: ', num2str(flux_id)]); end
        end
        
        % Step 5: Mass conservation check
        mass_conserve_check=Store_mat(i,1,t)-Store1(i)-Inflow(i)+Outflow_norm(i)+Outflow_evap(i);
        % If the deviation from mass conservation is over 0.0001%, warn
        if abs(mass_conserve_check)>Store_mat(i,1,t)/10^3 && abs(mass_conserve_check)>mass_check_value*1000
            disp(['WARNING: mass conserve error at time step: ',num2str(t), ' store: ', num2str(i),' error: ',num2str(mass_conserve_check)]); end
                
        % Step 6: Model the isotope concentration
        %   Check storage conditions to determine how it is modeled
        %       - Is there stored water?
        %       - Is there evaporation?
        %       - Is there transient throughflow?

        % Check for stored water:
        if Store_mat(i,1,t)>mass_check_value
            % Storage is wet
            % Calculate the 'simple mixed' isoStore and concentration
            isomixonly_conc=((isoStore1(i)+isoInflow(i))/(1+(Outflow_norm(i)+Outflow_evap(i))/Store_mat(i,1,t)))/Store_mat(i,1,t);
         
            % Check for evaporation and get the evaporation concentration if needed
            if Outflow_evap(i)>0
                E_ET_split=E_ET_split/Outflow_evap(i); %Weighted E/ET split
                % Use simple mix as source concentration:
                [ConcE,delstar]=craig_gordon_2H(isomixonly_conc, atm_mat(Store_mat(i,2,t),1:3,t), E_ET_split);
                % Find time 2 isotope store and concentration with fractionating evaporation:
                isoStore2(i)=(isoStore1(i)+isoInflow(i)-ConcE*Outflow_evap(i))/(1+Outflow_norm(i)/Store_mat(i,1,t));
                isoStore_conc(i,t+1)=isoStore2(i)/Store_mat(i,1,t);
                % Check fractionation limit if there is evaporation:
                %   Evaporation should not enrich beyond del*
                if isoStore_conc(i,t+1)>delta2conc_2H(delstar)
                    if isomixonly_conc>delta2conc_2H(delstar)
                        % if the source water is already more enriched than delstar
                        % there should be no further fractionation:
                        %   the final concentration is the 'simple mixed' value
                        isoStore_conc(i,t+1)=isomixonly_conc;
                        isoStore2(i)=isoStore_conc(i,t+1)*Store_mat(i,1,t);
                    else
                        % Evaporation can enrich, but only up to the
                        % limit based on local atm limits
                        isoStore_conc(i,t+1)=delta2conc_2H(delstar);
                        isoStore2(i)=isoStore_conc(i,t+1)*Store_mat(i,1,t);
                    end
                end
            else
                %No evaporation and 'simple mixed' is final concentration
                isoStore_conc(i,t+1)=isomixonly_conc;
                isoStore2(i)=isoStore_conc(i,t+1)*Store_mat(i,1,t);
            end        
            
        elseif Inflow(i)>0 || Store1(i)>0
            % No water left in storage but still either old water leaving or new water passed through
            isoStore_conc(i,t+1)=(isoInflow(i)+isoStore1(i))/(Inflow(i)+Store1(i));
            isoStore2(i)=isoStore_conc(i,t+1)*Store_mat(i,1,t);
        else
            % Storage was dry, is dry and stayed dry throughout dt
            isoStore_conc(i,t+1)=0; %no water, no concentration: -1000 per mil
            isoStore2(i)=isoStore_conc(i,t+1)*Store_mat(i,1,t);
        end
    end
    
    % Step 7: Update the storages for the next step
    Store1=Store_mat(:,1,t);
    isoStore1=isoStore2;
end

% Final return: all storages in delta format
isoStore_delta=conc2delta_2H(isoStore_conc);
end

function [ConcE,delstar] = craig_gordon_2H(source_conc,atm,evapFrac)
%This function returns the concentration of heavy isotopes in evaporated
%water vapor from a specified source, estimated based on atmospheric data
% Constants here are for 2H
% Equation references:  Horita and Wesolowski (1994), Gat & Levy (1978) and Gat (1981)

% Values for various factors listed here, if changes wanted or needed:
%     Scaled turbulence:
%     nexp=0.5 turbulent; =2/3 laminar; =1 static (soils/leaves)
%     fee=0.88 for large lake like the Great Lakes (Gat et al, 1994)
%     fee=1 for small water bodies where evap does not perturn ambient moisture
%     Cd=25.115 per mil (Gibson, Pisa)

% Scaling turbulence factor for MA simulations:
% if the weighted E/ET fraction is less than 50%, assume static soil water
% with higher E/ET, increase turbulence factor up to turbulent value for 100% E

nexp=1-(evapFrac-0.5);
if nexp>1; nexp=1; end
fee=1;
Cd=25.115;

% Convert atmosphere inputs to correct internal units
Temp_K=atm(1)+273.15; %assume input deg C, need deg K
delp=atm(3); %assume delta per mil, delta per mil needed
RH=atm(2)/100; %assume input % rel H, need fractional RH

% liquid to vapour equilibrium fractionation factor
% Horita and Wesolowski equation (1994)
alphastar=exp((1158.8*(Temp_K^3/10^9)-1620.1*(Temp_K^2/10^6)+794.84*(Temp_K/10^3)-161.04+2.9992*(10^9/Temp_K^3))/1000.);
estar=alphastar-1;

% Atmospheric vapour (delta a) is in equilibrium with isotopes in precipitation
% divide by alphastar since large scales (alphastar.ne.1)

dela=(delp/1000-estar)/alphastar;

% Source water (i.e. the water that is evaporating) in delta form (delta l)
dell=conc2delta_2H(source_conc); %per mil?

% Kinetic fractionation:
ekin=nexp*Cd*fee*(1-RH)/1000;

% Desication composition (delta *) and evaporation composition (delta e) 
% Gat & Levy (1978) and Gat (1981)
delstar=(RH*dela+ekin+estar/alphastar)/(RH-ekin-estar/alphastar)*1000;
dele=((dell-estar*1000)/alphastar-RH*dela*1000-ekin*1000)/(1-RH+ekin);

%Fractionation doesn't happen at very high RH, in such cases set the
%   evaporation composition to the source water composition
if dele>dell; dele=dell; end

ConcE=delta2conc_2H(dele);

end

function [delta] = conc2delta_2H(conc)
%Convert delta per mille to concentration for 2H
Rstd=0.00015576; %For internal accuracy using V-SMOW
R=conc./(1-conc);
delta=(R./Rstd-1)*1000;
end

function [conc] = delta2conc_2H(delta)
%Convert concentration to delta per mille for 2H
Rstd=0.00015576; %For internal accuracy using V-SMOW
R=(delta./1000+1)*Rstd;
conc = R./(R+1);
end