% A very basic hydrologic model for isotope simulation testing
%Three storages, 6 fluxes
%
%Init:
S_surface(1)=2;
S_ground(1)=10;
S_river(1)=5;
%Precip input
Precip=[0 5 3 0 0 12 5 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

for t=1:length(Precip)
        %Surface
    %Add precipitation
    S_surface(t+1)=S_surface(t)+Precip(t);
    %Conditional recharge
    if S_surface(t+1)>5
        Recharge(t)=5;
    else
        Recharge(t)=S_surface(t+1);
    end
    S_surface(t+1)=S_surface(t+1)-Recharge(t);
    %Conditional runoff
    if S_surface(t+1)>10
        Runoff(t)=(S_surface(t+1)-10)*0.5;
    else
        Runoff(t)=0;
    end
    S_surface(t+1)=S_surface(t+1)-Runoff(t);
    
    %Ground
    S_ground(t+1)=S_ground(t)+Recharge(t);
    Baseflow(t)=S_ground(t+1)*0.02;
    S_ground(t+1)=S_ground(t+1)-Baseflow(t);
    %River
    S_river(t+1)=S_river(t)+Runoff(t)+Baseflow(t);
    %Evaporation conditional (more E when no rain, lower E when low)
    if S_river(t+1)>1
        if Precip(t)>0
            Evap(t)=0.1;
        else
            Evap(t)=1;
        end
    else
        Evap(t)=S_river(t+1)*0.9;
    end
    S_river(t+1)=S_river(t+1)-Evap(t);
    %Outflow conditional
    if S_river(t+1)>0.2
        Outflow(t)=0.1*S_river(t+1)^1.5;
        if Outflow(t)>S_river(t+1)-10&&Outflow(t)>50
            Outflow(t)=S_river(t+1)-10; %when river is very high flow is still less than storage
        end
    else
        Outflow(t)=0; %no flow at low levels
    end
    S_river(t+1)=S_river(t+1)-Outflow(t);
end
    
    
