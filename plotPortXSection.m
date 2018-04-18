function plotPortXSection( port,ind )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

switch port.type
    case "circ"
        D=port.FinalDiameter;
        fw=port.fuelweb(ind);
        
        %plot the Outer Boundary
        pos = [-D/2 -D/2 D D];
        rectangle('Position',pos,'Curvature',[1 1],'EdgeColor',[0,0,0],'LineWidth',3)
        hold on
        
        %color the fuel area
        pos = [-D/2 -D/2 D D];
        rectangle('Position',pos,'Curvature',[1 1],'FaceColor',[0.9,0.9,0.9])
        
        
        %color the combusion chamber area
        pos = [-(D/2-fw) -(D/2-fw) D-2*fw D-2*fw];
        rectangle('Position',pos,'Curvature',[1 1],'FaceColor',[1,1,1])
        
        axis equal
        grid on
        title('Cross Section of Fuel Grain [meters]')
        
        
    case "Dport"
        D_outer = port.D_outer;
        fw = port.fuelweb(ind);
        tau = port.tau;
        
        %plot outer wall
        pos = [-D_outer/2 -D_outer/2 D_outer D_outer];
        rectangle('Position',pos,'Curvature',[1 1],'EdgeColor',[0,0,0],'LineWidth',3)
        hold on
        
        
        %plot fuel grain
        pos = [-D_outer/2 -D_outer/2 D_outer D_outer];
        rectangle('Position',pos,'Curvature',[1 1],'FaceColor',[0.9,0.9,0.9])
        
        %plot combustion area
        a=D_outer/2-fw;
        thetas = linspace(-acos((fw+tau)/a),acos((fw+tau)/a))';
        x=[a*cos(thetas) -a*cos(thetas)];
        y=[a*sin(thetas) a*sin(thetas)];
        
        patch(x,y,[1,1,1])
        
        
        %plot central bar (without caring about the curvature
        pos = [-tau -D_outer/2 2*tau D_outer];
        rectangle('Position',pos,'FaceColor',[0,0,0],'EdgeColor',[0,0,0],'LineWidth',3)
        
        axis equal
        axis on
        grid on
        title('Cross Section of Fuel Grain [meters]')
        
end

end

