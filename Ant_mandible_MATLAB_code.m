% % Calculate an estimate for the resultant force exerted on the apodeme from the muscle fibres.
    clear all
    Landmarks = readtable('atta30_landmarks.csv');
    
%Data below assume zero change in x displacement
    L = 4.46764*10^(-4);

%calculation of the opening angle when mandibles in a 'closed' state
    Mandible_Tooth = [(Landmarks.X(9)-Landmarks.X(11));(Landmarks.Y(9)-Landmarks.Y(11));(Landmarks.Z(9)-Landmarks.Z(11))]
    Head_spike_RL =  [(Landmarks.X(2)-Landmarks.X(1));(Landmarks.Y(2)-Landmarks.Y(1));(Landmarks.Z(2)-Landmarks.Z(1))]

    Unit_Mandible_Tooth = (1/norm(Mandible_Tooth))*Mandible_Tooth;
    Unit_Head_spike_RL = (1/norm(Head_spike_RL))*Head_spike_RL;

%closed mandible angle at 51 degrees, varies from 51-111 degrees.
    Beta = round(180-(acosd(dot(Unit_Head_spike_RL,Unit_Mandible_Tooth)/norm(Unit_Head_spike_RL))))

%% load data of closed mandible scan into workplace

    Data_fibre_length_orientation = readtable('atta30_cm_f_l_dat.csv');
    Orientation_vector = readtable('atta30_apo_rot.csv');
    
%direction of apodeme
    u = [0.04571567; -0.2369074; 0.97045606]; 
   
%adjust for voxel resolution of data and calculation of original pennation angle for all muscle fibers
  for i= 1:976
    Data_fibre_length_orientation.Length(i) = (2.44*10^(-6))*Data_fibre_length_orientation.Length(i);
    a = [ Data_fibre_length_orientation.X_dir(i);Data_fibre_length_orientation.Y_dir(i); Data_fibre_length_orientation.Z_dir(i)];
    Data_fibre_length_orientation.org_pennation(i) = acosd(dot(a,u)/norm(a));
  end

%% Calculate the muscle fibre length and pennation at each angle

% Do Example data between 51-111 degrees where beta = 51 degrees.
for i = 1:1:60
    New_fibre_data.Degrees(i) = (i + Beta);
end

% Calculate muscle fibre lengths and pennation angles at over 51-111 degrees opening angles with 1 degree intervals using equation 2

for j=1:976;
    for i = 1:1:60
        
         penn_org = Data_fibre_length_orientation.org_pennation(j);  
         FL = Data_fibre_length_orientation.Length(j);

         New_fibre_data.Pennation_angle(j,i) = atand((sind(penn_org) * FL)/((cosd(penn_org) * FL )+( sind(New_fibre_data.Degrees(i) - Beta) * L)));

         New_fibre_data.New_Length(j,i) = (sind(penn_org)/sind(New_fibre_data.Pennation_angle(j,i))) * FL

    end
end

%% 
for k=1:1:4
    for j=1:976;
        %set original fibre length
        FL = New_fibre_data.New_Length(j,((10*(k-1))+1));

        %assuming that the cross section area of fibre constant with contraction, calculate maximum force exerted by an individual muscle fibre
        Av_fibre_diameter = 27.7920369361307*10^(-6);
        Av_area = pi * ((Av_fibre_diameter/2)^2);

        Max_musc_stress = 300*(10^(3));
        Max_musc_force = Max_musc_stress *Av_area;

        %graph of force distrubution relative to the fibres stretch
        x = [0.2:0.01:1.8];
        y = Max_musc_force*((normpdf(x,1,0.2))/(normpdf(1,1,0.2)));
        figure(1)
        plot(x,y);


            for i =1:60
                %calculate indivudal fiber muscle force for every opening angle in direction of fibre
                Force_data.total_force(j,i) = Max_musc_force*((normpdf((New_fibre_data.New_Length(j,i)/FL),1,0.2))/(normpdf(1,1,0.2)));
            end
    end

    %% Calculate average muscle fibre length and pennation angle for all muscle fibers


    for i=1:1:60
        mean_pen = 0;
        total_pen = 0;
        mean_length = 0;
        total_length = 0;

        for j=1:1:976
            total_pen = total_pen + New_fibre_data.Pennation_angle(j,i);
            total_length = total_length + New_fibre_data.New_Length(j,i);
        end

        mean_pen = total_pen/976;
        mean_length = total_length/976;

        New_fibre_data.mean_pennation(1,i) = mean_pen;
        New_fibre_data.mean_length(1,i) = mean_length;
    end

%examples of fibre length and pennation
    x = [51:1:110]

    figure(2);
    plot(x,New_fibre_data.Pennation_angle(1,:))
    hold on
    plot(x,New_fibre_data.Pennation_angle(2,:))
    hold on
    plot(x,New_fibre_data.Pennation_angle(3,:))
    hold on
    plot(x,New_fibre_data.Pennation_angle(4,:))
    xlabel('Opening angle (degrees)')
    ylabel('Pennation angle (degrees)')
    legend('Example 1','Example 2','Example 3','Example 4')

    figure(3);
    curve_fx = polyfit(x,New_fibre_data.mean_pennation(1,:),3);
    v = polyval(curve_fx, x);
    plot(x,v)
    hold on
    xlabel('Opening angle (degrees)')
    ylabel('Mean pennation angle (degrees)')

    figure(4);
    plot(x,New_fibre_data.mean_length(1,:))
    hold on
    xlabel('Opening angle (degrees)')
    ylabel('Mean fibre length (m)')


    figure(5);
    plot(x,New_fibre_data.New_Length(1,:))
    hold on
    plot(x,New_fibre_data.New_Length(2,:))
    hold on
    plot(x,New_fibre_data.New_Length(3,:))
    hold on
    plot(x,New_fibre_data.New_Length(4,:))
    xlabel('Opening angle (degrees)')
    ylabel('Muscle fibre length (m)')
    legend('Example 1','Example 2','Example 3','Example 4')

%% Calculate the resultant force of each fibre in the direction of the apodeme at each opening angle & summing them together for the total force

for i=1:60;
    force_val = 0;
    for j=1:976;
        force_val = force_val +(Force_data.total_force(j,i)*cosd(New_fibre_data.Pennation_angle(j,i)))
        Force_data.total_musc_force.F_apo_tot(k,i) = force_val;
    end
end


end

%plot of opening angle vs force on apodeme with different muscle relaxation states
    figure(6)
    x = [51:1:110]
    plot(x,Force_data.total_musc_force.F_apo_tot)
    xlabel('Opening angle (degrees)')
    ylabel('Fibre force (N)')
    lgd = legend('51^{\circ}','61^{\circ}','71^{\circ}','81^{\circ}');
    title(lgd,'Opening angle where the muscle fibre assumed to be in a relaxed state','FontSize',12)

%% Import all data for experimental and simulated, and calculate fitted plot

    X = [70 75 80 85 90 95 100]

%Bite force data for muscle fibers in relaxed state at 51 degrees opening angle
    Y_51 = (1*10^(3))*[0.026 0.015 0.008 0.003 0.0013 0.00065 0.0003]
    curve_Y51 = polyfit(X,Y_51,3);
    Y51 = polyval(curve_Y51, X);

%Bite force data for muscle fibers in relaxed state at 61 degrees opening angle
    Y_61 = (1*10^(3))*[0.07 0.05 0.04 0.0225 0.013 0.0085 0.006]
    curve_Y61 = polyfit(X,Y_61,3);
    Y61 = polyval(curve_Y61, X);
    
%Bite force data for muscle fibers in relaxed state at 71 degrees opening angle
    Y_71 = (1*10^(3))*[0.09 0.093 0.085 0.062 0.05 0.035 0.022]
    curve_Y71 = polyfit(X,Y_71,3);
    Y71 = polyval(curve_Y71, X);
    
%Bite force data for muscle fibers in relaxed state at 81 degrees opening angle
    Y_81 = (1*10^(3))*[0.072 0.09 0.105 0.092 0.085 0.08 0.05]
    curve_Y81 = polyfit(X,Y_81,3);
    Y81 = polyval(curve_Y81, X);

lgd = legend('51^{\circ}','61^{\circ}','71^{\circ}','81^{\circ}');
title(lgd,'Opening angle where the muscle fibre assumed to be in a relaxed state','FontSize',12)

%import experimental data
    FY = [127.7567874 107.3009942 93.76202146 107.4646294 97.19460767 87.6435746 62.54694011 80.96293756 70.55299195 71.82906882 70.70277438 49.84752392 99.69389743 90.49635801 42.64806523 74.65072496 84.80984824 70.88070402 42.87143646 59.62945578 56.93956541 46.75476142 43.9330685 89.79058248 72.06557752 45.19480485 40.88102837 33.84326668 41.50074519 41.06746323 22.26206078] 
    FX = [68.65628042 74.24648614 74.83540757 75.74408635 76.02880226 76.71068344 77.40583484 79.72706221 80.89688831 81.93227467 82.09711637 82.22164399 82.452114 82.60517958 83.0868093 83.77364348 85.64581885 85.68592824 86.45444147 87.07272641 89.07632592 90.40750133 90.90073545 91.41950818 91.83460639 92.74913449 95.42965423 98.33216909 99.36377348 99.67138333 100.9221504]

    
%% create data and line for comparison of normalised  bite force - opening angle relationship (simulated & experimental data)

Comp_y51 = Y_51/max(Y_51);
Comp_y61 = Y_61/max(Y_61);
Comp_y71 = Y_71/max(Y_71);
Comp_y81 = Y_51/max(Y_81);


Comp_fy = FY/max(FY);
curve_fx = polyfit(FX,Comp_fy,3);
v = polyval(curve_fx, FX);
curve_ffx = polyfit(FX,FY,3);
P = polyval(curve_ffx, FX);



curve_x51 = polyfit(X,Comp_y51,3);
v51 = polyval(curve_x51, X);
curve_x61 = polyfit(X,Comp_y61,3);
v61 = polyval(curve_x61, X);
curve_x71 = polyfit(X,Comp_y71,3);
v71 = polyval(curve_x71, X);
curve_x81 = polyfit(X,Comp_y81,3);
v81 = polyval(curve_x81, X);


%% Create plots
figure(1)
tiledlayout(2,1);

%Plot describing the bite force - opening angle relationship (simulated data)
nexttile
    plot(X,Y51);
    hold on;
    plot(X,Y61);
    hold on
    plot(X,Y71);
    hold on
    plot(X,Y81);
    hold on
    title('Bite force - opening angle (Simulink Model)')
    xlabel('Opening angle (degrees)')
    ylabel('Modelled Bite force (mN)')
    lgd = legend('51^{\circ}','61^{\circ}','71^{\circ}','81^{\circ}');
    title(lgd,'Opening angle where the muscle fibre assumed to be in a relaxed state','FontSize',12)
    
%Plot describing the bite force - opening angle relationship (experimental data)
nexttile
    plot(FX,FY,'x')
    hold on
    plot(FX,P)
    title('Bite force - opening angle (Experimental)')
    xlabel('Opening angle (degrees)')
    ylabel('Modelled Bite force (mN)')
    legend('Experimental data','Experimental data - line of best fit')
hold off

%Comparison of normalised  bite force - opening angle relationship (simulated & experimental data)

figure(2)
    plot(X,v51)
    hold on
    plot(X,v61)
    hold on
    plot(X,v71)
    hold on
    plot(X,v81)
    hold on
    plot(FX,v)
    hold on

    title('Comparison of normalised bite force - opening angle (Simulated vs Experimental)')
    xlabel('Opening angle (degrees)')
    ylabel('Normalised modelled Bite force (N)')
    legend('Simulation data-unstretched muscle at 51^{\circ}','Simulation data-unstretched muscle at 61^{\circ}','Simulation data-unstretched muscle at 71^{\circ}','Simulation data-unstretched muscle at 81^{\circ}','Experimental data - line of best fit')

    
 