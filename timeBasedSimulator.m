function [rocketPerformance, rocketSim] = timeBasedSimulator(universalConstants,rocketDesign,regRateParams, deltaT,thresh)


%% extract variables

%universal constants
g0    = universalConstants.g0;
bar   = universalConstants.bar;
P_amb = universalConstants.P_amb;

a = regRateParams.a;
n = regRateParams.n;
m = regRateParams.m;

etac = rocketDesign.etac;
lambda = rocketDesign.lambda;

Lp = rocketDesign.Lp;
rho_fuel = rocketDesign.rho_fuel;

P_cc=rocketDesign.P_cc;

mdot_oxinit = rocketDesign.mdot_oxinit;
mdot_ox = mdot_oxinit; %assumed constant throughout burn
m_ox = rocketDesign.m_ox;

A_exit = rocketDesign.A_exit;
A_throat = rocketDesign.A_throat;
expansionRatio = rocketDesign.expansionRatio;

rocketSim.port=rocketDesign.port;

if rocketSim.port.type == "circ"
    rocketSim.port.diameter=rocketSim.port.IntialDiameter;
    rocketSim.port.fuelweb = rocketSim.port.InitialFuelWeb;
end

if rocketSim.port.type == "Dport"
    rocketSim.port.fuelweb = rocketSim.port.InitialFuelWeb;
end



%% peform simulation

qburnfin = 0; %parameter that determines whether the burn is completed or not (0 when not completed, 1 when completed)

t=0;
ti=1; %time index


while qburnfin == 0


    % 2. determine fuel regression rate & % 3. determine total mass flow rate
    
    % Calculating area, perimeter of current port diameter
    [A_port, P_port] = portCalculator(rocketSim.port);
    
    G_ox(ti) = mdot_ox/(A_port);
    
    G_prop(ti) = G_ox(ti); %as a starting point
    
    val(1) = G_prop(ti); %stored temporarily to check for convergence
    
    for loopind = 1:1000
        
        S(ti) = P_port*Lp;
        
        rdot(ti) = a*G_prop(ti)^n*Lp^m;
        
        
        mdot_fuel(ti) = rdot(ti)*rho_fuel*S(ti);
        
        
        G_fuel(ti) = mdot_fuel(ti)./(A_port); %#ok<*SAGROW>
        
        G_prop(ti) = G_ox(ti)+G_fuel(ti);
        
        val(loopind+1) = G_prop(ti);
        
        if abs(val(loopind)-val(loopind+1))/val(loopind)<thresh
            break;
        end
        
    end
    
    mdot_prop(ti)=mdot_ox+mdot_fuel(ti); %[kg/s]
    
    % 4. determine thermochemistry
    clear val
    
    val(1) = P_cc(max([1,ti-1]));
    
    for loopind = 1:1000
        
        P_cc(ti) = val(loopind);
        
        OF(ti)=mdot_ox/mdot_fuel(ti);
        
        [T_flame(ti), gamma(ti), m_mol(ti), R(ti), c_star(ti)] = thermochem(OF(ti),P_cc(ti),etac); % runs a script called thermochem to determine these numbers.
        
        T_stag(ti) = T_flame(ti); %Approximation! Need better formula for stagnation temp of flow
        
        
        % 5. determine eng perf
        
        P_cc(ti) = mdot_prop(ti)*c_star(ti)/A_throat; %assume there is no throat erosion for now - should be updated
        
        val(loopind+1) = P_cc(ti);
        
        if abs(val(loopind)-val(loopind+1))/val(loopind)<thresh
            break;
        end
        
    end
    
    P_stag(ti) = P_cc(ti); %approximation! need to check if better formula exists.
    
    [M_exit(ti), Tratio(ti), Pratio(ti), rhoratio(ti), ~] = flowisentropic(gamma(ti), expansionRatio, 'sup');
    
    T_exit(ti) = T_stag(ti)*Tratio(ti);
    P_exit(ti) = P_stag(ti)*Pratio(ti);
    
    v_exit(ti) = M_exit(ti)*sqrt(gamma(ti)*R(ti)*T_exit(ti));
    
    F(ti) = lambda*(mdot_prop(ti)*v_exit(ti) + (P_exit(ti)-P_amb)*A_exit); %thrust
    
    Isp(ti) = F(ti)/(mdot_prop(ti)*g0);
    
    % 6. update geometry
    switch rocketSim.port.type
        case "circ"
            rocketSim.port.fuelweb(ti+1)=rocketSim.port.fuelweb(ti)-rdot(ti)*deltaT;
            rocketSim.port.diameter(ti+1)=rocketSim.port.diameter(ti)+rdot(ti)*deltaT;
        case "Dport"
            rocketSim.port.fuelweb(ti+1)=rocketSim.port.fuelweb(ti)-rdot(ti)*deltaT;
    end
    
    %throat erosion
    
           
    ti=ti+1;        %increment the time index
    t=t+deltaT;     %increment the time
    
    % 8. check if burn is completed, and go to next step
    
    %qburnfin = checkburnfin(); % a full function could be written to
    %perform this check properly.At the moment it just checks for burn
    %time completion.
    
    %Things to check:
    % (1) chamber pressure < ox tank pressure
    % (2) Choked flow
    % (3) fuel is not entirely eaten up
    
    if mdot_ox*t > m_ox %oxidiser used up
        qburnfin=1;
        disp('Ox used up');
    else
        if rocketSim.port.fuelweb(end)<=0 %using only condition 3 for now
            qburnfin=1;
            disp('Fuel Web used up')
        else
            qburnfin=0;
        end
    end
        
    
    
    
end


%% store results

rocketSim.F = F;
rocketSim.Isp= Isp;
rocketSim.mdot_fuel=mdot_fuel;
rocketSim.mdot_prop=mdot_prop;
rocketSim.c_star=c_star;
rocketSim.G_ox=G_ox;
rocketSim.G_fuel=G_fuel;
rocketSim.G_prop= G_prop;
rocketSim.M_exit=M_exit;
rocketSim.OF = OF;
rocketSim.P_cc = P_cc;
rocketSim.P_exit=P_exit;


rocketPerformance.I_total_result = sum(F)*deltaT;

rocketPerformance.t_burn_result = t;

rocketPerformance.F_init_result  = F(1);

rocketPerformance.F_avg_result  = mean(F);

rocketPerformance.Isp_avg = mean(Isp);

rocketPerformance.m_f_total = sum(mdot_fuel)*deltaT;

rocketPerformance.m_ox_total = mdot_ox*rocketPerformance.t_burn_result;

rocketPerformance.m_prop_total = sum(mdot_prop)*deltaT;
















end

