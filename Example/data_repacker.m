function [Store_init,Store_mat,Flux_mat,atm_mat,time1,time2] = data_repacker(S_ground,S_river,S_surface,Baseflow,Evap,Outflow,Precip,Recharge,Runoff)
%Repackage data for complete model run of the VBH for MAITsim
%   Create init matrix, time vectors
%   Repack storages and fluxes
% Store init:
Store_init(1,1) = S_surface(1);
Store_init(2,1) = S_ground(1);
Store_init(3,1) = S_river(1);
% isoStore init (delta):
Store_init(1,2) = -12;
Store_init(2,2) = -20;
Store_init(3,2) = -15;

% Matrix packing loop
for n=2:length(S_river)
    %Time:
    time1(n-1)=n-1;
    time2(n-1)=n;

    %Storage:
    %store2
    Store_mat (1,1,n-1)= S_surface(n);
    Store_mat (2,1,n-1)= S_ground(n);
    Store_mat (3,1,n-1)= S_river(n);
    %atmosphere
    Store_mat (1,2,n-1)= 1;
    Store_mat (2,2,n-1)= 1;
    Store_mat (3,2,n-1)= 1;
    
    %Atmosphere (only one for this model):
    atm_mat(1,1,n-1)= 15; %Temperature
    atm_mat(1,2,n-1)= 88; %RH
    atm_mat(1,3,n-1)= -10; %dP
    
    %Fluxes:
    % flux, source, destination, evap frac
    Flux_mat(1,1:4,n-1)= [Precip(n-1) 0 1 0];
    Flux_mat(2,1:4,n-1)= [Evap(n-1) 3 0 1];
    Flux_mat(3,1:4,n-1)= [Runoff(n-1) 1 3 0];
    Flux_mat(4,1:4,n-1)= [Recharge(n-1) 1 2 0];
    Flux_mat(5,1:4,n-1)= [Baseflow(n-1) 2 3 0];
    Flux_mat(6,1:4,n-1)= [Outflow(n-1) 3 0 0];
    
end
end

