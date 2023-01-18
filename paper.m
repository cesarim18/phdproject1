function paper
%% LOAD PARAMETERS

%% Strain Vs Stress
    clf
    N_Ro = 2000;
    Ro = linspace(0.8,1.1,N_Ro);
    Ro = reshape(Ro,N_Ro,1);

    N_Strain = 2000;
    Strain = linspace(1.05,1.25,N_Strain);
    Strain = reshape(Strain,N_Strain,1);
    
    % con grosor = 1.25 mm y Go = 40 kPa, se ajustan los datos
    
    
    grosor = 2.5*0.5/10.0;% 0.5 mm
    Go = 20.0; % 10 kPa
    Po = 9.0; % 3 kPa
    beta = 2.0; % adimensional
    P = Go*(Strain.^beta - 1);
    % Go = (4.0/3.0)*EY*hd/rd;
    for k=1:200:N_Ro
        RoK = Ro(k,1);
        R = RoK.*Strain;
        R_I = R - 0.5*grosor;
        R_E = R + 0.5*grosor;
        Stress_Go = 2.0*P./((R./R_I).^2 - (R./R_E).^2);
        Stress_Po = 2.0*Po./((R./R_I).^2 - (R./R_E).^2);
        plot(Strain(1:35:N_Strain),Stress_Go(1:35:N_Strain) + Stress_Po(1:35:N_Strain),'--')
        hold on
        pause(0.01)
        k
    end
    Data = load('Extracted.csv');
    plot(Data(:,1),Data(:,2)./10000.0,'o')
         xlabel('Strain','FontSize',24)
         hold on
         ylabel('Stress (kPa)','FontSize',24)
         hold on

    print('Strain-Stress_Extracted_Data','-depsc',figure(1))
stop

%% Strain Vs Stress
    clf
    N_Ro = 100;
    Ro = linspace(0.8,1.1,N_Ro);
    Ro = reshape(Ro,N_Ro,1);

    N_Strain = 2000;
    Strain = linspace(1.05,1.4,N_Strain);
    Strain = reshape(Strain,N_Strain,1);
    
    % con grosor = 1.25 mm y Go = 40 kPa, se ajustan los datos
    
    
    grosor = 2.5*0.5/10.0;% 0.5 mm
    Go = 20.0; % 10 kPa
    Po = 7.5; % 3 kPa
    beta = 2.0; % adimensional
    P = Go*(Strain.^beta - 1);
    
    for k=1:10:N_Ro
        RoK = Ro(k,1);
        R = RoK.*Strain;
        R_I = R - 0.5*grosor;
        R_E = R + 0.5*grosor;
        Stress_Go = 2.0*P./((R./R_I).^2 - (R./R_E).^2);
        Stress_Po = 2.0*Po./((R./R_I).^2 - (R./R_E).^2);
        plot(Strain(1:20:N_Strain),Stress_Go(1:20:N_Strain) + Stress_Po(1:20:N_Strain),'.')
        hold on
        pause(0.01)
        k
    end
    Data = load('Activation_Data.csv');
    plot(Data(:,1),100.0*Data(:,2),'o')
    hold on

    print('Strain-Stress_Activation_Data','-depsc',figure(1))
    
stop    
%% valores que aproximan los datos
%    clf
%    figure(1)
    Data = load('Control_Data.csv')
    Strain = linspace(0.8,1.5,2000);
    h = 0.5; % milimetros
    R0 = 0.8; % centimetros
    beta = 2.0; % adimensional
    G0 = 10.0; % kilopascales
    P0 = 3.0; %kilopascales
    for i=1:2000
        RRI(1,i) = Strain(1,i)/(Strain(1,i) - 0.5*h/R0/10.0);
        RRE(1,i) = Strain(1,i)/(Strain(1,i) + 0.5*h/R0/10.0);
        Stress_P0(1,i) = 2*P0/(RRI(1,i)^2 - RRE(1,i)^2);
        Pressure(1,i) = G0*(Strain(1,i)^beta - 1);
        Stress_G0(1,i) = 2*Pressure(1,i)/(RRI(1,i)^2 - RRE(1,i)^2);
    end
    %plot(Strain(1,:),Stress_P0(1,:),'--');
    %hold on
    %plot(Strain(1,:),Stress_G0(1,:),'--');
    %hold on
    plot(Strain(1,:),Stress_P0(1,:) + Stress_G0(1,:),'--');
    hold on
   plot(Data(:,1),10.0*Data(:,2),'.')
   hold on
%    axis([1.09 1.23 50 200])
%    hold on
    xlabel('Strain','FontSize',24)
    hold on
    ylabel('Stress','FontSize',24)
    hold on
    print('Strain-Stress_Approach','-depsc',figure(1))

        
    stop
%%


    N_s = 200+4;
    N_th = 3*60;
	delta_s = ((7.0357+0.8+0.9+6.4737+15.2+1.8+0.7+0.7+4.3+4.3)/100.0)/(N_s-4)
	delta_th = 2.0*pi/N_th

	Vessel_x = zeros(N_s-4,N_th+1);
	Vessel_y = zeros(N_s-4,N_th+1);
	Vessel_z = zeros(N_s-4,N_th+1);
	Vessel_r = zeros(N_s-4,N_th+1);

	y = zeros(N_s,1);
	L = (7.0357+0.8+0.9+6.4737+15.2+1.8+0.7+0.7+4.3+4.3)/100.0;
    y0 = 0;
    for j = 1:N_s
		j_shift=j-2;
		y_j=y0+delta_s*(j_shift-2);
		y(j,1) = y_j;
    end
    
    y0 = delta_s*(-2.0);
    axis_L = linspace(y0,L-y0,N_s);
    axis_P = linspace(0,2*pi,N_th+1);
	fr = 2;
	Tfinal = 1.0;
	movies = 101;
	Ite = movies-1;
	Tmov = Tfinal/Ite;

    for i = 1:movies
        axis_T(i) = (i-1)*Tmov;
    end
    
    
%    PT = load('Example/Ex203/Files/PT.dat');
    P = load('Example/Ex203/Files/P.dat');
    R = load('Example/Ex203/Files/R.dat');
    R_0 = load('Example/Ex203/IC/IC_Ro.dat');
    G_0 = load('Example/Ex203/IC/IC_Go.dat');
    P0 = 3.0; %kilopascales
    grosor = 2.5*0.5/10.0;% 1.25 mm
    for j = 1:size(R,1)
        R_I(j,1) = 100.0*R(j,1) - 0.5*grosor;
        R_E(j,1) = 100.0*R(j,1) + 0.5*grosor;
    end
    mov_R = size(R,1)/(N_s)/(N_th+1)
    mov_P = size(P,1)/(N_s)/(N_th+1)
    
    
    for i = 1:mov_R
        for j = 1:N_s*(N_th+1)
            k = j+(i-1)*N_s*(N_th+1);
            Strain_2D(k,1) = R(k,1)/R_0(j,1);
            Stress_2D(k,1) = 2*(1.05*P(k,1) + P0)/((100.0*R(k,1)/R_I(k,1))^2 - (100.0*R(k,1)/R_E(k,1))^2);
        end
    end
 	R = reshape(R,N_s,N_th+1,mov_R);
 	P = reshape(P,N_s,N_th+1,mov_R);
%    P(:,1,1)
%    stop
    R_I = reshape(R_I,N_s,N_th+1,mov_R);
	R_E = reshape(R_E,N_s,N_th+1,mov_R);
    G_0 = reshape(G_0,N_s,N_th+1,2);
    R_0 = reshape(R_0,N_s,N_th+1,2);

    min_S = min(Stress_2D(:,1));
    max_S = max(Stress_2D(:,1));
    Strain_2D = reshape(Strain_2D,N_s,N_th+1,mov_R);
    Stress_2D = reshape(Stress_2D,N_s,N_th+1,mov_R);
    Stress_2D(:,1,1)
    %stop
    Data = load('Extracted.csv')

    angle1 = 1;
    angle2 = 16;
    angle3 = 31;
    angle4 = 46;

%     subplot(231)
%         title(['Begin the second CC'])
%         hold on
%         mov = 1;
%         plot(Strain_2D(3:1:N_s-2,angle1,mov),Stress_2D(3:1:N_s-2,angle1,mov),'.')
%         hold on
%         plot(Strain_2D(3:1:N_s-2,angle2,mov),Stress_2D(3:1:N_s-2,angle2,mov),'.')
%         hold on
%         mov = 21;
%         plot(Strain_2D(3:1:N_s-2,angle1,mov),Stress_2D(3:1:N_s-2,angle1,mov),'.')
%         hold on
%         plot(Strain_2D(3:1:N_s-2,angle2,mov),Stress_2D(3:1:N_s-2,angle2,mov),'.')
%         hold on
%         %plot(Strain_2D(3:1:N_s-2,angle3,mov),Stress_2D(3:1:N_s-2,angle3,mov),'.')
%         %hold on
%         %plot(Strain_2D(3:1:N_s-2,angle4,mov),Stress_2D(3:1:N_s-2,angle4,mov),'.')
%         %hold on
%         plot(Data(:,1),Data(:,2)./10000.0,'.')
%         hold on
%         axis([0.95 1.25 min_S 200])
%         hold on
%         xlabel('Strain','FontSize',24)
%         hold on
%         ylabel('Stress','FontSize',24)
%         hold on
%     
%     subplot(232)
%         title(['Begin the third CC'])
%         hold on
%         mov = 201;
%         plot(Strain_2D(3:1:N_s-2,angle1,mov),Stress_2D(3:1:N_s-2,angle1,mov),'.')
%         hold on
%         plot(Strain_2D(3:1:N_s-2,angle2,mov),Stress_2D(3:1:N_s-2,angle2,mov),'.')
%         hold on
%         mov = 221;
%         plot(Strain_2D(3:1:N_s-2,angle1,mov),Stress_2D(3:1:N_s-2,angle1,mov),'.')
%         hold on
%         plot(Strain_2D(3:1:N_s-2,angle2,mov),Stress_2D(3:1:N_s-2,angle2,mov),'.')
%         hold on
%         %plot(Strain_2D(3:1:N_s-2,angle3,mov),Stress_2D(3:1:N_s-2,angle3,mov),'.')
%         %hold on
%         %plot(Strain_2D(3:1:N_s-2,angle4,mov),Stress_2D(3:1:N_s-2,angle4,mov),'.')
%         %hold on
%         plot(Data(:,1),Data(:,2)./10000.0,'.')
%         hold on
%         axis([0.95 1.25 min_S 200])
%         hold on
%         xlabel('Strain','FontSize',24)
%         hold on
%         ylabel('Stress','FontSize',24)
%         hold on
% 
%     subplot(233)
%         title(['Begin the fourth CC'])
%         hold on
%         mov = 301;
%         plot(Strain_2D(3:1:N_s-2,angle1,mov),Stress_2D(3:1:N_s-2,angle1,mov),'.')
%         hold on
%         plot(Strain_2D(3:1:N_s-2,angle2,mov),Stress_2D(3:1:N_s-2,angle2,mov),'.')
%         hold on
%         mov = 321;
%         plot(Strain_2D(3:1:N_s-2,angle1,mov),Stress_2D(3:1:N_s-2,angle1,mov),'.')
%         hold on
%         plot(Strain_2D(3:1:N_s-2,angle2,mov),Stress_2D(3:1:N_s-2,angle2,mov),'.')
%         hold on
%         %plot(Strain_2D(3:1:N_s-2,angle3,mov),Stress_2D(3:1:N_s-2,angle3,mov),'.')
%         %hold on
%         %plot(Strain_2D(3:1:N_s-2,angle4,mov),Stress_2D(3:1:N_s-2,angle4,mov),'.')
%         %hold on
%         plot(Data(:,1),Data(:,2)./10000.0,'.')
%         hold on
%         axis([0.95 1.25 min_S 200])
%         hold on
%         xlabel('Strain','FontSize',24)
%         hold on
%         ylabel('Stress','FontSize',24)
%         hold on
% 
%         print('Example/Ex203/Videos/Strain-Stress','-depsc',figure(1))
%stop    

    Movie = VideoWriter('Example/Ex203/Videos/Stress-Strain_Horiz.avi');
	Movie.Quality = 100;
	Movie.FrameRate = 3;
    for mov=1:mov_R
        hold off
%        plot(Strain(3:1:N_s-2,angle,mov),Stress(3:1:N_s-2,angle,mov))
%        hold on
         plot(Strain_2D(3:1:N_s-2,angle1,mov),Stress_2D(3:1:N_s-2,angle1,mov),'.')
         hold on
%          plot(Strain_2D(3:1:N_s-2,angle1,mov+100),Stress_2D(3:1:N_s-2,angle1,mov+100),'.')
%          hold on
%          plot(Strain_2D(3:1:N_s-2,angle1,mov+200),Stress_2D(3:1:N_s-2,angle1,mov+200),'.')
%          hold on
%          plot(Strain_2D(3:1:N_s-2,angle1,mov+300),Stress_2D(3:1:N_s-2,angle1,mov+300),'.')
%          hold on
%          plot(Strain_2D(3:1:N_s-2,angle1,mov+400),Stress_2D(3:1:N_s-2,angle1,mov+400),'.')
%          hold on
%        plot(Strain_2D(3:1:N_s-2,angle2,mov),Stress_2D(3:1:N_s-2,angle2,mov),'.')
%        hold on
%        plot(Strain_2D(3:1:N_s-2,angle3,mov),Stress_2D(3:1:N_s-2,angle3,mov),'.')
%        hold on
%        plot(Strain_2D(3:1:N_s-2,angle4,mov),Stress_2D(3:1:N_s-2,angle4,mov),'.')
%        hold on
         plot(Strain_2D(3:1:N_s-2,angle2,mov),Stress_2D(3:1:N_s-2,angle2,mov),'.')
         hold on
%          plot(Strain_2D(3:1:N_s-2,angle2,mov+100),Stress_2D(3:1:N_s-2,angle2,mov+100),'.')
%          hold on
%          plot(Strain_2D(3:1:N_s-2,angle2,mov+200),Stress_2D(3:1:N_s-2,angle2,mov+200),'.')
%          hold on
%          plot(Strain_2D(3:1:N_s-2,angle2,mov+300),Stress_2D(3:1:N_s-2,angle2,mov+300),'.')
%          hold on
%          plot(Strain_2D(3:1:N_s-2,angle2,mov+400),Stress_2D(3:1:N_s-2,angle2,mov+400),'.')
%          hold on
        %plot(Strain_2D(3:1:N_s-2,angle2,mov+100),Stress_2D(3:1:N_s-2,angle2,mov+100),'.')
        %hold on
        %plot(Strain_2D(3:1:N_s-2,angle2,mov+200),Stress_2D(3:1:N_s-2,angle2,mov+200),'.')
        %hold on
        plot(Data(:,1),Data(:,2)./10000.0,'.')
        hold on
        axis([0.95 1.25 min_S max_S])
        hold on
        xlabel('Strain','FontSize',24)
        hold on
        ylabel('Stress','FontSize',24)
        hold on
        time = (mov-1)*Tmov;
        %legend('0-1 second 0°','1-2 second 0°','2-3 second 0°','0-1 second 90°','1-2 second 90°','2-3 second 90°','Data Collected')
        %legend('0°','90°','180°','270°','Data Collected')
        title(['Strain vs. Stress at t = ', num2str(time),' seconds for'])
        mov
	  	open(Movie)
	 	F = getframe(figure(1));
	 	writeVideo(Movie,F);
        pause(0.1);
    end
    close(Movie);
%stop


    U = load('Example/Ex203/Files/u.dat');
    Vx = load('Example/Ex203/Files/Vx.dat');
    Vy = load('Example/Ex203/Files/Vy.dat');
    Vz = load('Example/Ex203/Files/Vz.dat');
	U = reshape(U,N_s,N_th+1,mov_R);
    Vx = reshape(Vx,N_s,N_th+1,mov_R);
	Vy = reshape(Vy,N_s,N_th+1,mov_R);
	Vz = reshape(Vz,N_s,N_th+1,mov_R);

    Vessel_x = zeros(N_s-4,N_th+1);
	Vessel_y = zeros(N_s-4,N_th+1);
	Vessel_z = zeros(N_s-4,N_th+1);
	Vessel_r = zeros(N_s-4,N_th+1);

    MovAvi1 = VideoWriter('Example/Ex203/Videos/3D_Strain.avi');
	MovAvi1.Quality = 100;
	MovAvi1.FrameRate = fr;
    Qv_s = 4;
    Qv_th = 3;
    for mov = 1:mov_R
        for k = 1:N_th+1
		    for j = 1:N_s-4
		        Vessel_x(j,k) = Vx(j+2,k,mov);
		        Vessel_y(j,k) = Vy(j+2,k,mov);
		        Vessel_z(j,k) = Vz(j+2,k,mov);
		        Vessel_r(j,k) = Strain_2D(j+2,k,mov);
		    end
		end
		clf
		figure(1) %splot = 
		COLOR = Vessel_r;
		surf(Vessel_x,Vessel_y,Vessel_z,COLOR); 
		c = jet(100);%jet(100);
		colormap(c);
		colorbar;
		caxis([1 1.15])
		%caxis([-1 1])
		%caxis([CC_UsL CC_UsR])
		%caxis([R2L R2R])
        shading interp;
		axis equal
		camlight('left')
		hold on
		%h1 = quiver3(Vx(3:Qv_s:end,1:Qv_th:end,mov),Vy(3:Qv_s:end,1:Qv_th:end,mov),Vz(3:Qv_s:end,1:Qv_th:end,mov),Velx(3:Qv_s:end,1:Qv_th:end,mov),Vely(3:Qv_s:end,1:Qv_th:end,mov),Velz(3:Qv_s:end,1:Qv_th:end,mov),'r');
        %axis([-0.02 0.1 -0.02 0.02 -0.32 0.07])
        %hold on
        xlabel('x in m','FontSize',12)
		%ylabel('y in m','FontSize',12)
		zlabel('z in m','FontSize',12)
		set(gca,'FontSize',12)
		Tmov1 = (mov-1)*Tmov;
		%title(['Radius Percentage, T = ', num2str(Tmov1),' seconds'])
		title(['Strain, T = ', num2str(Tmov1),' seconds'])
		%title(['Tangential Velocity, T = ', num2str(Tmov1),' seconds'])
		%view(-90,0)
		mov
	  	open(MovAvi1)
	 	F = getframe(figure(1));
	 	writeVideo(MovAvi1,F);
		%pause(0.1)
    end
    close(MovAvi1);
stop
    MovAvi0 = VideoWriter('Example/Ex203/Videos/1D_Us.avi');
    MovAvi0.Quality = 100;
    MovAvi0.FrameRate = fr;
    for mov = 1:101
        %mov = 21
        clf
        figure(1)
        Tmov1 = (mov-1)*Tmov;
        xlabel('s in m');
        ylabel('Axial velocity in ms^-1');
        title(['Axial Velocity U, T = ', num2str(Tmov1),' seconds']);
        hold on
        plot(axis_L,U(:,3*N_th/4 + 1,mov),'r');
        hold on
        plot(axis_L,U(:,N_th/4 + 1,mov),'b');
        hold on
        %axis([y0 L-y0 CC_UsL CC_UsR])
        axis([axis_L(1) axis_L(N_s) -1 1])
        hold on
        mov
        %pause
 	  	open(MovAvi0)
 	 	F = getframe(figure(1));
 	 	writeVideo(MovAvi0,F);
    end
    close(MovAvi0);
    
    
    
    
 stop
    for i=1:mov_R
        clf
        subplot(231)
        plot(Strain_2D(3:1:N_s-2,1,i))
        hold on
        subplot(232)
        plot(Stress_2D(3:1:N_s-2,1,i))
        hold on
        subplot(233)
        plot(Strain_2D(3:1:N_s-2,1,i),Stress_2D(3:1:N_s-2,1,i))
        hold on
        pause(0.1)
    end
    stop    
%% -- 25/05/22 19:01:47 --%
    Data = load('Extracted.csv')
    b1=Data(:,2)\Data(:,1)
    yCalc1 = b1*Data(:,1)
    scatter(Data(:,1),Data(:,2)/10000.0)
    hold on
    %plot(Data(:,1),0.00001*yCalc1)
    X = [ones(length(Data(:,1)),1) Data(:,1)]
    b = X\Data(:,2)
    yCalc2 = X*b
    plot(Data(:,1),yCalc2/10000.0,'--')


%% contraste con 2D modelo  
end
