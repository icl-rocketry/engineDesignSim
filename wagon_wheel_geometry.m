function [ L_wagon,r_total,r_central ] = wagon_wheel_geometry( n_max,A_port,mdot_fuelinit,GO_init,a,m,n,rho_fuel,m_f )
% calculates the fuel grain geometry for a wagon wheel shape. note
%that with an aerospike, the central radius denotes the radius of the
%central support column of the spike

%last edited: Will H. feb 16 2018

%initialisation
L_wagon = zeros(1,n_max-2);
web_wagon = zeros(1,n_max-2);
r_central = zeros(1,n_max-2);
r_total = zeros(1,n_max-2);

nport = 3;

spike_support = 0.01; %[m] thickness of walls required to support the aerospike central column

%Loops from minimum number of ports n=3 up to the specified maximum number
while(nport <= n_max)
    
    num_ports = nport;
    
    if nport < 7 %under 7, curvature of the outer wall of each port is significant, and thus the geometry changes
        
        sidelength = sqrt((A_port*nport)/pi);
        
        perimeter_init = 2*sidelength + (pi*2*sidelength)/num_ports;
        
        theta_port = pi/num_ports;
        
    else %above 7 ports, the ports are effectively triangular, and the geometry can be approximated as such
       
        theta_port = pi/num_ports;

        height_init = sqrt(A_port/(tan(theta_port)));

        base_init = 2*height_init*tan(theta_port);

        sidelength = height_init/cos(theta_port);

        perimeter_init = 2*sidelength + base_init;
        
    end
    
    GF_wagon_init = mdot_fuelinit/(num_ports*A_port);

    L_wagon(nport-2) = ((mdot_fuelinit)/(num_ports*rho_fuel*a*((GO_init + GF_wagon_init)^n)*perimeter_init))^(1/(m+1));

    web_poly = [1,perimeter_init/pi, -(m_f)/(num_ports*rho_fuel*L_wagon(nport-2)*pi)];

    web_roots = roots(web_poly);

    for i = 1:2
        
        if web_roots(i) > 0
            
            web_wagon = web_roots(i);
            
        end
        
    end
    
    r_central(nport-2) = (web_wagon/sin(theta_port)) - web_wagon;
    
    if nport < 7
        r_total(nport-2) = sidelength + 2*(web_wagon+spike_support) + r_central(nport-2);
    else
        r_total(nport-2) = height_init + 2*(web_wagon+spike_support) + r_central(nport-2);
    end
    
    nport = nport+1;
    
end



end

