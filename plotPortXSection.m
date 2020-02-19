function plotPortXSection( port,ind )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

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

end

