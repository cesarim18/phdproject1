function blood
	N_s = 100+4;
    N_th = 60;
    R_IC = load('Example/Ex201/IC/IC_Ro.dat');
    G_IC = load('Example/Ex201/IC/IC_Go.dat');
    R_IC = reshape(R_IC,N_s,N_th+1,2)
    G_IC = reshape(G_IC,N_s,N_th+1,2)

    clf
    figure(1)
    subplot(221)
    plot(100.*axis_L,100.*R_IC(:,1,1),'--r','Linewidth',4)
    hold on
    plot(100.*axis_L,100.*R_IC(:,N_th/4,1),'b','Linewidth',4)
    hold on
    axis([0 100*L 0.8 1.7])
    set(gca,'FontSize',24)
    hold on
    xlabel('Arc length (cm)','Fontsize',24)
    ylabel('R_o (cm)','Fontsize',24)
    legend('0 rad','\pi/2 rad')
    
    subplot(222)
    plot(100.*axis_L,1.05.*G_IC(:,1,1),'k','Linewidth',4)
    hold on
    axis([0 100*L 47 68])
    hold on
    set(gca,'FontSize',24)
    xlabel('Arc length (cm)','Fontsize',24)
    ylabel('G_o (kPa)','Fontsize',24)

    subplot(223)
    plot(axis_L,100.*R_IC(:,1,1),'r')
    hold on
    plot(axis_L,-100.*R_IC(:,1,1),'b')
    hold on
    axis([0 L -1.7 1.7])
    hold on

    subplot(224)
    plot(axis_L,0.5*(100.*R_IC(:,1,1)).^2,'k')
    hold on
%% LOAD GEOMETRY
    CC_Vx = load('Example/Ex201/Files/Vx.dat');
	mov_CC_Vx = size(CC_Vx,1)/(N_s)/(N_th+1)
	CC_Vx = reshape(CC_Vx,N_s,(N_th+1),mov_CC_Vx);

    CC_Vy = load('Example/Ex201/Files/Vy.dat');
	mov_CC_Vy = size(CC_Vy,1)/(N_s)/(N_th+1)
	CC_Vy = reshape(CC_Vy,N_s,(N_th+1),mov_CC_Vy);

	CC_Vz = load('Example/Ex201/Files/Vz.dat');
	mov_CC_Vz = size(CC_Vz,1)/(N_s)/(N_th+1)
	CC_Vz = reshape(CC_Vz,N_s,(N_th+1),mov_CC_Vz);

%% LOAD VELOCITY FIELD
    CC_Velx = load('Example/Ex201/Files/Velx.dat');
	mov_CC_Velx = size(CC_Velx,1)/(N_s)/(N_th+1)
	CC_Velx = reshape(CC_Velx,N_s,(N_th+1),mov_CC_Velx);

    CC_Vely = load('Example/Ex201/Files/Vely.dat');
	mov_CC_Vely = size(CC_Vely,1)/(N_s)/(N_th+1)
	CC_Vely = reshape(CC_Vely,N_s,(N_th+1),mov_CC_Vely);

    CC_Velz = load('Example/Ex201/Files/Velz.dat');
	mov_CC_Velz = size(CC_Velz,1)/(N_s)/(N_th+1)
	CC_Velz = reshape(CC_Velz,N_s,(N_th+1),mov_CC_Velz);

%% AXIAL VELOCITY
    CC_Us = load('Example/Ex98/Files/u.dat');
	CC_UsL = min(CC_Us(:,1))
	CC_UsR = max(CC_Us(:,1))
	if(CC_UsR - CC_UsL == 0.0)
		CC_UsL = -10^-15
		CC_UsR = -CC_UsL
    end
 	mov_CC_Us = size(CC_Us,1)/(N_s)/(N_th+1)
	CC_Us = reshape(CC_Us,N_s,N_th+1,mov_CC_Us);

    Ex2_Comp = -CC_Us + Ex2_Us;
    Ex3_Comp = -CC_Us + Ex3_Us;
    Ex4_Comp = -CC_Us + Ex4_Us;
    Ex5_Comp = -CC_Us + Ex5_Us;
    

    clf
    figure(1)
    subplot(221)
    j = N_s/6;
    s_pos = 100*(j/N_s)*L;
    k = N_th/4+1;
    for i = 1:mov_CC_Us
        Pjk(i) = CC_Us(j,k,i);
    end
    title(['1/6'])
    hold on
    plot(axis_T,Pjk,'k','Linewidth',3)
    k = 3*N_th/4+1;
    for i = 1:mov_CC_Us
        Pjk3(i) = CC_Us(j,k,i);
    end
    plot(axis_T,Pjk3,'-.r','Linewidth',3)
    axis([0 2 -1 1])
    hold on
    legend('90°','270°')

    subplot(222)
    j = N_s/2;
    s_pos = 100*(j/N_s)*L;
    k = N_th/4+1;
    for i = 1:mov_CC_Us
        Pjk(i) = CC_Us(j,k,i);
    end
    title(['1/2'])
    hold on
    plot(axis_T,Pjk,'k','Linewidth',3)
    k = 3*N_th/4+1;
    for i = 1:mov_CC_Us
        Pjk3(i) = CC_Us(j,k,i);
    end
    plot(axis_T,Pjk3,'-.r','Linewidth',3)
    axis([0 2 -1 1])
    hold on
    legend('90°','270°')


%% RADIUS
    CC_R = load('Example/Ex201/Files/R.dat');
	CC_RL = min(CC_R(:,1))
	CC_RR = max(CC_R(:,1))
	mov_CC_R = size(CC_R,1)/(N_s)/(N_th+1)
	CC_R = reshape(CC_R,N_s,N_th+1,mov_CC_R);
    size(R_IC)/N_s/(N_th+1)
    for i = 1:mov_CC_R
        for j = 1:N_s*(N_th+1)
            CC_R2(j+(i-1)*N_s*(N_th+1),1) = CC_R(j+(i-1)*N_s*(N_th+1),1)/R_IC(j,1);
        end
    end    
    %R_IC = reshape(R_IC,N_s,N_th+1,1);    
    CC_R2L = min(CC_R2(:,1))
	CC_R2R = max(CC_R2(:,1))
	CC_R2 = reshape(CC_R2,N_s,N_th+1,mov_CC_R);
    R_IC = reshape(R_IC,N_s,N_th+1,1);

%% Pressure
    CC_P = load('Example/Ex201/Files/P.dat');
	CC_PL = min(CC_P(:,1))
	CC_PR = max(CC_P(:,1))
	mov_CC_P = size(CC_P,1)/(N_s)/(N_th+1)
    CC_P = reshape(CC_P,N_s,N_th+1,mov_CC_P);
    Go_IC = load('Example/Ex201/IC/IC_Go.dat');
    Go_IC = reshape(Go_IC,N_s,N_th+1,1);

%% TANGENTIAL VELOCITY
    CC_Ut = load('Example/Ex201/Files/Ut.dat');
    CC_UtL = min(CC_Ut(:,1))
	CC_UtR = max(CC_Ut(:,1))
  	mov_CC_Ut = size(CC_Ut,1)/(N_s)/(N_th+1)
	CC_Ut = reshape(CC_Ut,N_s,N_th+1,mov_CC_Ut);


%% COMPARATION WITH CARDIAC CYCLE

    MovAvi0 = VideoWriter('Example/Ex98/Videos/1D_Us.avi');
    MovAvi0.Quality = 100;
    MovAvi0.FrameRate = fr;
    for mov = 1:mov_CC_Us
        %mov = 21
        clf
        figure(1)
        Tmov1 = (mov-1)*Tmov;
        xlabel('s in m');
        ylabel('Axial velocity in ms^-1');
        title(['Axial Velocity U_s, T = ', num2str(Tmov1),' seconds']);
        hold on
        plot(axis_L,CC_Us(:,3*N_th/4 + 1,mov),'r');
        hold on
        plot(axis_L,CC_Us(:,N_th/4 + 1,mov),'b');
        hold on
        %axis([y0 L-y0 CC_UsL CC_UsR])
        axis([axis_L(3) axis_L(N_s-2) CC_UsL CC_UsR])
        hold on
        mov
        %pause
 	  	open(MovAvi0)
 	 	F = getframe(figure(1));
 	 	writeVideo(MovAvi0,F);
    end
    close(MovAvi0);

    MovAvi1 = VideoWriter('Example/Ex98/Videos/3D_CC_Us.avi');
	MovAvi1.Quality = 100;
	MovAvi1.FrameRate = fr;
    Qv_s = 4;
    Qv_th = 3;
    for mov = 1:mov_CC_Vx
    %mov=mov_R
        for k = 1:N_th+1
		    for j = 1:N_s-4
		        Vessel_x(j,k) = CC_Vx(j+2,k,mov);
		        Vessel_y(j,k) = CC_Vy(j+2,k,mov);
		        Vessel_z(j,k) = CC_Vz(j+2,k,mov);
		        Vessel_r(j,k) = CC_Us(j+2,k,mov);
		    end
		end
		clf
		figure(1) %splot = 
		COLOR = Vessel_r;
		surf(Vessel_x,Vessel_y,Vessel_z,COLOR); 
		c = jet(100);%jet(100);
		colormap(c);
		colorbar;
		caxis([-0.5 1])
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
		title(['Axial Velocity, T = ', num2str(Tmov1),' seconds'])
		%title(['Tangential Velocity, T = ', num2str(Tmov1),' seconds'])
		view(-90,0)
		mov
	  	open(MovAvi1)
	 	F = getframe(figure(1));
	 	writeVideo(MovAvi1,F);
		%pause(0.1)
    end
    close(MovAvi1);

    clf
    figure(1)
    j = N_s/2;
    s_pos = 100*(j/N_s)*L;
    k = N_th/4+1;
    for i = 1:mov_CC_Us
        Pjk(i) = CC_Us(j,k,i);
    end
    subplot(221)
    title(['Angle 90'])
    hold on
    plot(axis_T,Pjk)

    k = N_th/2+1;
    for i = 1:mov_CC_Us
        Pjk2(i) = CC_Us(j,k,i);
    end
    subplot(222)
    title(['Angle 180'])
    hold on
    plot(axis_T,Pjk2)

    k = 3*N_th/4+1;
    for i = 1:mov_CC_Us
        Pjk3(i) = CC_Us(j,k,i);
    end
    subplot(223)
    title(['Angle 270'])
    hold on
    plot(axis_T,Pjk3)

    k = N_th+1;
    for i = 1:mov_CC_Us
        Pjk4(i) = CC_Us(j,k,i);
    end
    subplot(224)
    title(['Angle 360'])
    hold on
    plot(axis_T,Pjk4)
    sgtitle(['Axial Velocity Profiles in time, s = ', num2str(s_pos),' cm'])
    hold on
    print('Example/Ex201/Videos/Waves_AxialVelocity','-depsc',figure(1))

    clf
    figure(1)
    j = N_s/2;
    s_pos = 100*(j/N_s)*L;
    k1 = N_th/4;
    k2 = N_th/2;
    k3 = 3*N_th/4;
    k4 = N_th;
    t1 = (k1/N_th)*360;
    t2 = (k2/N_th)*360;
    t3 = (k3/N_th)*360;
    t4 = (k4/N_th)*360;
    for i = 1:mov_CC_Us
        Us_k1(i) = CC_Us(j,k1,i);
        Us_k2(i) = CC_Us(j,k2,i);
        Us_k3(i) = CC_Us(j,k3,i);
        Us_k4(i) = CC_Us(j,k4,i);
    end
    plot(axis_T,Us_k1,'b--','Linewidth',4)
    hold on
    plot(axis_T,Us_k2,'r','Linewidth',4)
    hold on
    plot(axis_T,Us_k3,'y','Linewidth',4)
    hold on
    plot(axis_T,Us_k4,'g--','Linewidth',4)
    hold on
    plot(axis_T,0.0.*Us_k4,'k','Linewidth',4)
    hold on
    axis([0 2 -0.2 0.8])
    hold on
    xlabel('Time (s)','FontSize',24)
    ylabel('Axial Velocity (m/s)','FontSize',24)
    title(['Profiles at s = ', num2str(s_pos),' cm'])
	set(gca,'FontSize',24)
    legend('90°','180°','270°','360°')
    print('Example/Ex201/Videos/Waves_AxialVelocity_','-depsc',figure(1))



    MovAvi1 = VideoWriter('Example/Ex201/Videos/3D_R.avi');
	MovAvi1.Quality = 100;
	MovAvi1.FrameRate = 2;
    Qv_s = 4;
    Qv_th = 3;
    for mov = 1:mov_CC_Vx
       % mov=1
        for k = 1:N_th+1
		    for j = 1:N_s
		        jk = j + N_s*(N_th+1-k);
		        Vessel_x(j,k) = CC_Vx(j,k,mov);
		        Vessel_y(j,k) = CC_Vy(j,k,mov);
		        Vessel_z(j,k) = CC_Vz(j,k,mov);
		        Vessel_r(j,k) = R(j,k,mov);
		    end
		end
		clf
		figure(1) %splot = 
		COLOR = Vessel_r;
		surf(Vessel_x,Vessel_y,Vessel_z,COLOR); 
		c = jet(100);%jet(100);
		colormap(c);
		colorbar;
		%caxis([UsL UsR])
		%caxis([CC_RL CC_RR])
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
		title(['Radius, T = ', num2str(Tmov1),' seconds'])
		%title(['Axial Velocity, T = ', num2str(Tmov1),' seconds'])
		%title(['Tangential Velocity, T = ', num2str(Tmov1),' seconds'])
		view(-118,16)
		mov
	  	open(MovAvi1)
	 	F = getframe(figure(1));
	 	writeVideo(MovAvi1,F);
		%pause(0.1)
    end
    close(MovAvi1);
    
    MovAvi0 = VideoWriter('Example/Ex95/Videos/1D_R.avi');
    MovAvi0.Quality = 100;
    MovAvi0.FrameRate = fr;
    for mov = 1:mov_R
%        mov = 1
        clf
        figure(1)
        Tmov1 = (mov-1)*Tmov;
        xlabel('s in m');
        ylabel('Radius');
        title(['Radius, T = ', num2str(Tmov1),' seconds']);
        hold on
        plot(axis_L,R(:,3*N_th/4,mov),'b');
        hold on
        plot(axis_L,-R(:,N_th/4,mov),'b');
        hold on
        plot(axis_L,R_I(:,3*N_th/4,mov),'g');
        hold on
        plot(axis_L,-R_I(:,N_th/4,mov),'g');
        hold on
        plot(axis_L,R_0(:,3*N_th/4,1),'r');
        hold on
        plot(axis_L,-R_0(:,N_th/4,1),'r');
        hold on
        %axis([y0 L-y0 -CC_RR CC_RR])
        %hold on
        mov
        %pause
 	  	open(MovAvi0)
 	 	F = getframe(figure(1));
 	 	writeVideo(MovAvi0,F);
    end
    close(MovAvi0);

    MovAvi0 = VideoWriter('Example/Ex201/Videos/1D_R2.avi');
    MovAvi0.Quality = 100;
    MovAvi0.FrameRate = fr;
    for mov = 1:mov_CC_R
        clf
        figure(1)
        Tmov1 = (mov-1)*Tmov;
        xlabel('s in m');
        ylabel('Radius');
        title(['Radius, T = ', num2str(Tmov1),' seconds']);
        hold on
        plot(axis_L,CC_R2(:,3*N_th/4,mov),'b');
        hold on
        plot(axis_L,-CC_R2(:,N_th/4,mov),'b');
        hold on
         plot(axis_L,R_IC(:,3*N_th/4,1)./R_IC(:,3*N_th/4,1),'r');
         hold on
         plot(axis_L,-R_IC(:,N_th/4,1)./R_IC(:,N_th/4,1),'r');
         hold on
        axis([y0 L-y0 -CC_R2R CC_R2R])
        hold on
        mov
        %pause
 	  	open(MovAvi0)
 	 	F = getframe(figure(1));
 	 	writeVideo(MovAvi0,F);
    end
    close(MovAvi0);

    clf
    figure(1)
    j = N_s/2;
    s_pos = 100*(j/N_s)*L;
    k1 = N_th/4;
    k2 = N_th/2;
    k3 = 3*N_th/4;
    k4 = N_th;
    for i = 1:mov_CC_R
        R_k1(i) = CC_R(j,k1,i);
        R_k2(i) = CC_R(j,k2,i);
        R_k3(i) = CC_R(j,k3,i);
        R_k4(i) = CC_R(j,k4,i);
    end
    plot(axis_T,R_k1,'b','Linewidth',4)
    hold on
    plot(axis_T,R_k2,'r','Linewidth',4)
    hold on
    plot(axis_T,R_k3,'y--','Linewidth',4)
    hold on
    plot(axis_T,R_k4,'g--','Linewidth',4)
    hold on
    plot(axis_T,0.0.*R_k4,'k','Linewidth',4)
    hold on
    %axis([0 2 -0.2 0.8])
    %hold on
    xlabel('Time (s)','FontSize',24)
    ylabel('Radius (m/s)','FontSize',24)
    title(['Profiles at s = ', num2str(s_pos),' cm'])
	set(gca,'FontSize',24)
    legend('90°','180°','270°','360°')
    print('Example/Ex201/Videos/Waves_Radius_','-depsc',figure(1))

    
    MovAvi1 = VideoWriter('Example/Ex201/Videos/3D_P.avi');
	MovAvi1.Quality = 100;
	MovAvi1.FrameRate = fr;
    Qv_s = 4;
    Qv_th = 3;
    for mov = 1:mov_CC_P
		for k = 1:N_th+1
		    k1 = mod(k-2+(N_th/2),N_th+1)+1;
		    for j = 1:N_s
		        jk = j + N_s*(N_th+1-k);
		        Vessel_x(j,k) = CC_Vx(j,k,mov);
		        Vessel_y(j,k) = CC_Vy(j,k,mov);
		        Vessel_z(j,k) = CC_Vz(j,k,mov);
		        Vessel_r(j,k) = CC_P(j,k,mov);
		    end
		end
		clf
		figure(1) %splot = 
		COLOR = Vessel_r;
		surf(Vessel_x,Vessel_y,Vessel_z,COLOR); 
		c = jet(100);%jet(100);
		colormap(c);
		colorbar;
		caxis([CC_PL CC_PR])
        shading interp;
		axis equal
		camlight('left')
		hold on
		%h1 = quiver3(Vx(3:Qv_s:end,1:Qv_th:end,mov),Vy(3:Qv_s:end,1:Qv_th:end,mov),Vz(3:Qv_s:end,1:Qv_th:end,mov),Velx(3:Qv_s:end,1:Qv_th:end,mov),Vely(3:Qv_s:end,1:Qv_th:end,mov),Velz(3:Qv_s:end,1:Qv_th:end,mov),'r');
        %axis([-0.02 0.1 -0.02 0.02 -0.32 0.07])
        %hold on
        xlabel('x in m','FontSize',12)
		zlabel('z in m','FontSize',12)
		set(gca,'FontSize',12)
		Tmov1 = (mov-1)*Tmov;
		title(['Pressure, T = ', num2str(Tmov1),' seconds'])
		view(-118,16)
		mov
	  	open(MovAvi1)
	 	F = getframe(figure(1));
	 	writeVideo(MovAvi1,F);
    end
    close(MovAvi1);

    MovAvi0 = VideoWriter('Example/Ex201/Videos/1D_P.avi');
    MovAvi0.Quality = 100;
    MovAvi0.FrameRate = 2;
    for mov = 1:mov_R
        clf
        figure(1)
        Tmov1 = (mov-1)*Tmov;
        xlabel('s in m');
        ylabel('Pressure');
        title(['Pressure, T = ', num2str(Tmov1),' seconds']);
        hold on
        plot(P(N_s/2,:,mov));
        hold on
        plot(P(N_s,:,mov));
        hold on
%        axis([axis_L(3) axis_L(N_s-2) CC_PL CC_PR])
%        hold on
        mov
        %pause
 	  	open(MovAvi0)
 	 	F = getframe(figure(1));
 	 	writeVideo(MovAvi0,F);
    end
    close(MovAvi0);

%% plots for P
    clf
    figure(1)
    j = N_s/4;
    s_pos = 100*(j/N_s)*L;
    k1 = N_th/4;
    k2 = N_th/2;
    k3 = 3*N_th/4;
    k4 = N_th;
    for i = 1:mov_CC_P
        P_k1(i) = CC_P(j,k1,i)*1.05;
        P_k2(i) = CC_P(j,k2,i)*1.05;
        P_k3(i) = CC_P(j,k3,i)*1.05;
        P_k4(i) = CC_P(j,k4,i)*1.05;
    end
    plot(axis_T,P_k1,'b','Linewidth',4)
    hold on
    plot(axis_T,P_k2,'r','Linewidth',4)
    hold on
    plot(axis_T,P_k3,'y--','Linewidth',4)
    hold on
    plot(axis_T,P_k4,'g--','Linewidth',4)
    hold on
    plot(axis_T,0.0.*P_k4,'k','Linewidth',4)
    hold on
    %axis([0 2 -0.2 0.8])
    %hold on
    xlabel('Time (s)','FontSize',24)
    ylabel('Pressure (kPa)','FontSize',24)
    title(['Profiles at s = ', num2str(s_pos),' cm'])
	set(gca,'FontSize',24)
    legend('90°','180°','270°','360°')
    print('Example/Ex201/Videos/Waves_Pressure','-depsc',figure(1))


    MovAvi1 = VideoWriter('Example/Ex201/Videos/3D_Ut.avi');
	MovAvi1.Quality = 100;
	MovAvi1.FrameRate = fr;
    Qv_s = 15;
    Qv_th = 5;
    for mov = 1:mov_CC_Ut
		for k = 1:N_th+1
		    k1 = mod(k-2+(N_th/2),N_th+1)+1;
		    for j = 1:N_s
		        jk = j + N_s*(N_th+1-k);
		        Vessel_x(j,k) = CC_Vx(j,k,mov);
		        Vessel_y(j,k) = CC_Vy(j,k,mov);
		        Vessel_z(j,k) = CC_Vz(j,k,mov);
		        Vessel_r(j,k) = CC_Ut(j,k,mov);
		    end
		end
		clf
		figure(1) %splot = 
		COLOR = Vessel_r;
		surf(Vessel_x,Vessel_y,Vessel_z,COLOR); 
		c = jet(100);%jet(100);
		colormap(c);
		colorbar;
		caxis([CC_UtL CC_UtR])
        shading interp;
		axis equal
		camlight('left')
		hold on
        xlabel('x in m','FontSize',12)
		zlabel('z in m','FontSize',12)
		set(gca,'FontSize',12)
		Tmov1 = (mov-1)*Tmov;
		title(['Tangential Velocity, T = ', num2str(Tmov1),' seconds'])
		view(-118,16)
		mov
	  	open(MovAvi1)
	 	F = getframe(figure(1));
	 	writeVideo(MovAvi1,F);
    end
    close(MovAvi1);

    
    clf
    figure(1)
    j = N_s/2;
    s_pos = 100*(j/N_s)*L;
    k1 = N_th/4;
    k2 = N_th/2;
    k3 = 3*N_th/4;
    k4 = N_th;
    for i = 1:mov_CC_Ut
        Ut_k1(i) = CC_Ut(j,k1,i);
        Ut_k2(i) = CC_Ut(j,k2,i);
        Ut_k3(i) = CC_Ut(j,k3,i);
        Ut_k4(i) = CC_Ut(j,k4,i);
    end
    plot(axis_T,Ut_k1,'b','Linewidth',4)
    hold on
    plot(axis_T,Ut_k2,'r','Linewidth',4)
    hold on
    plot(axis_T,Ut_k3,'y--','Linewidth',4)
    hold on
    plot(axis_T,Ut_k4,'g--','Linewidth',4)
    hold on
    plot(axis_T,0.0.*Ut_k4,'k','Linewidth',4)
    hold on
    %axis([0 2 -0.2 0.8])
    %hold on
    xlabel('Time (s)','FontSize',24)
    ylabel('Tangential Velocity (m/s)','FontSize',24)
    title(['Profiles at s = ', num2str(s_pos),' cm'])
	set(gca,'FontSize',24)
    legend('90°','180°','270°','360°')
    print('Example/Ex201/Videos/Waves_TangentialVelocity','-depsc',figure(1))


    
    
%% Initial Conditions
    clf
    figure(1)
    subplot(142)
    mov = 1
        for k = 1:N_th+1
		    for j = 1:N_s-4
		        Vessel_x(j,k) = 100.*CC_Vx(j+2,k,mov);
		        Vessel_y(j,k) = 100.*CC_Vy(j+2,k,mov);
		        Vessel_z(j,k) = 100.*CC_Vz(j+2,k,mov);
%                Vessel_r(j,k) = CC_R(j+2,k,mov);
		    end
		end
		surf(Vessel_x,Vessel_y,Vessel_z); 
		c = gray(100);%jet(100);
		colormap(c);
		%colorbar;
        shading interp;
		axis equal
		camlight('left')
		hold on
		%h1 = quiver3(Vx(3:Qv_s:end,1:Qv_th:end,mov),Vy(3:Qv_s:end,1:Qv_th:end,mov),Vz(3:Qv_s:end,1:Qv_th:end,mov),Velx(3:Qv_s:end,1:Qv_th:end,mov),Vely(3:Qv_s:end,1:Qv_th:end,mov),Velz(3:Qv_s:end,1:Qv_th:end,mov),'r');
        %axis([-0.02 0.1 -0.02 0.02 -0.32 0.07])
        %hold on
        xlabel('x (cm)','FontSize',24)
		ylabel('y (cm)','FontSize',24)
		zlabel('z (cm)','FontSize',24)
		set(gca,'FontSize',24)
		Tmov1 = (mov-1)*Tmov;
		title(['Vessel'])
		view(-60,30)

        subplot(132)
    mov = 1
        for k = 1:N_th+1
		    for j = 1:N_s-4
		        Vessel_x(j,k) = CC_Vx(j+2,k,mov);
		        Vessel_y(j,k) = CC_Vy(j+2,k,mov);
		        Vessel_z(j,k) = CC_Vz(j+2,k,mov);
%                Vessel_r(j,k) = CC_R(j+2,k,mov);
		    end
		end
		surf(Vessel_x,Vessel_y,Vessel_z); 
		c = gray(100);%jet(100);
		colormap(c);
		%colorbar;
        shading interp;
		axis equal
		camlight('left')
		hold on
		%h1 = quiver3(Vx(3:Qv_s:end,1:Qv_th:end,mov),Vy(3:Qv_s:end,1:Qv_th:end,mov),Vz(3:Qv_s:end,1:Qv_th:end,mov),Velx(3:Qv_s:end,1:Qv_th:end,mov),Vely(3:Qv_s:end,1:Qv_th:end,mov),Velz(3:Qv_s:end,1:Qv_th:end,mov),'r');
        %axis([-0.02 0.1 -0.02 0.02 -0.32 0.07])
        %hold on
        xlabel('x (m)','FontSize',24)
		ylabel('y (m)','FontSize',24)
		zlabel('z (m)','FontSize',24)
		set(gca,'FontSize',24)
		Tmov1 = (mov-1)*Tmov;
		title(['Vessel'])
		view(-118,16)


    clf
    figure(1)
    font_size=20
    subplot(324)
    k0 = 1;
    k1 = N_th/4+1;
    k2 = N_th/2+1;
    k3 = 3*N_th/4+1;
    k4 = N_th+1;
    plot(100.*axis_L,100.*R_IC(:,k0,1),'b','Linewidth',4)
    hold on
    plot(100.*axis_L,100.*R_IC(:,k1,1),'-r','Linewidth',4)
    hold on
    axis([100.*axis_L(3) 100.*axis_L(N_s-2) 0.8 1.7])
    %hold on
    %xlabel('Time (s)','FontSize',24)
    ylabel('R_o (cm)','FontSize',font_size)
    %title(['Profiles at s = ', num2str(s_pos),' cm'])
	set(gca,'FontSize',font_size)
    legend('0 rad','\pi/2 rad')

    subplot(326)
    plot(100.*axis_L,G_IC(:,k1,1)*1.05,'k','Linewidth',4)
    hold on
    axis([100.*axis_L(3) 100.*axis_L(N_s-2) 45 70])
    %hold on
    xlabel('Arc-length (cm)','FontSize',font_size)
    ylabel('G_o (kPa)','FontSize',font_size)
    %title(['Profiles at s = ', num2str(s_pos),' cm'])
	set(gca,'FontSize',font_size)

% 	fileu1=load('Example/Ex201/Parameters/dt.dat');
% 	maxdt = max(fileu1(:,1));
%     mindt = min(fileu1(:,1));
%     tamdt = length(fileu1);
%     tamdt
%     fileu5=load('Example/Ex201/Parameters/u_gh.dat');
% 	maxgh = max(fileu5(:,1));
%     mingh = min(fileu5(:,1));
% 	tamgh = length(fileu5)
% 
%     T = zeros(tamdt,1);% linspactwie(0,fileu1(tamdt:1),tamdt+1)
%     T(1,1) = fileu1(1,1);
%     %stop
%     for i=2:tamdt
%         T(i) = T(i-1) + fileu1(i,1);
%         %for j=1:i
%         %T(i) = T(i) + fileu1(j,1);
%         %end
%     end
%     Tfinal=T(tamdt,1)

    subplot(322)
	plot(T,fileu5,'k','Linewidth',4)
    hold on
	axis([0 1 -0.2 1.2])
    hold on
	xlabel('Time (s)');
	ylabel({'Velocity at left';'boundary (m/s)'})
	set(gca,'FontSize',font_size)

%%
    clf
    figure(1)

    subplot(224)
    j = N_s/2;
    s_pos = 100*(j/N_s)*L;
    k0 = 1;
    k1 = N_th/4+1;
    k2 = N_th/2+1;
    k3 = 3*N_th/4+1;
    k4 = N_th+1;
    for i = 1:mov_CC_Ut
        Ut_k0(i) = CC_Ut(j,k0,i);
        Ut_k1(i) = CC_Ut(j,k1,i);
        Ut_k2(i) = CC_Ut(j,k2,i);
        Ut_k3(i) = CC_Ut(j,k3,i);
        Ut_k4(i) = CC_Ut(j,k4,i);
    end
    plot(axis_T,Ut_k0,'b','Linewidth',4)
    hold on
    plot(axis_T,Ut_k1,'k','Linewidth',4)
    hold on
    plot(axis_T,Ut_k2,'r--','Linewidth',4)
    hold on
    plot(axis_T,Ut_k3,'g--','Linewidth',4)
    hold on
    %axis([0 2 -0.2 0.8])
    %hold on
    xlabel('Time (s)','FontSize',24)
    ylabel('Tangential Velocity (m/s)','FontSize',24)
    %title(['Profiles at s = ', num2str(s_pos),' cm'])
	set(gca,'FontSize',24)
    %legend('90°','180°','270°','360°')


    subplot(221)
    j = N_s/2;
    s_pos = 100*(j/N_s)*L;
    k0 = 1;
    k1 = N_th/4+1;
    k2 = N_th/2+1;
    k3 = 3*N_th/4+1;
    k4 = N_th+1;
    for i = 1:mov_CC_R
        R_k0(i) = CC_R(j,k0,i)*100;
        R_k1(i) = CC_R(j,k1,i)*100;
        R_k2(i) = CC_R(j,k2,i)*100;
        R_k3(i) = CC_R(j,k3,i)*100;
        R_k4(i) = CC_R(j,k4,i)*100;
    end
    plot(axis_T,R_k0,'b','Linewidth',4)
    hold on
    plot(axis_T,R_k1,'k','Linewidth',4)
    hold on
    plot(axis_T,R_k2,'r--','Linewidth',4)
    hold on
    plot(axis_T,R_k3,'g--','Linewidth',4)
    hold on
    %axis([0 2 -0.2 0.8])
    %hold on
    %xlabel('Time (s)','FontSize',24)
    ylabel('Radius (m)','FontSize',24)
    %title(['Profiles at s = ', num2str(s_pos),' cm'])
	set(gca,'FontSize',24)
    legend('0 rad','\pi/2 rad','\pi rad','3\pi/2 rad')

    subplot(222)
    j = N_s/2;
    s_pos = 100*(j/N_s)*L;
    k0 = 1;
    k1 = N_th/4+1;
    k2 = N_th/2+1;
    k3 = 3*N_th/4+1;
    k4 = N_th+1;
    t1 = (k1/N_th)*360;
    t2 = (k2/N_th)*360;
    t3 = (k3/N_th)*360;
    t4 = (k4/N_th)*360;
    for i = 1:mov_CC_Us
        Us_k0(i) = CC_Us(j,k0,i);
        Us_k1(i) = CC_Us(j,k1,i);
        Us_k2(i) = CC_Us(j,k2,i);
        Us_k3(i) = CC_Us(j,k3,i);
        Us_k4(i) = CC_Us(j,k4,i);
    end
    plot(axis_T,Us_k0,'b','Linewidth',4)
    hold on
    plot(axis_T,Us_k1,'k','Linewidth',4)
    hold on
    plot(axis_T,Us_k2,'r--','Linewidth',4)
    hold on
    plot(axis_T,Us_k3,'g--','Linewidth',4)
    hold on
    axis([0 2 -0.2 0.8])
    hold on
    %xlabel('Time (s)','FontSize',24)
    ylabel('Axial Velocity (m/s)','FontSize',24)
    %title(['Profiles at s = ', num2str(s_pos),' cm'])
	set(gca,'FontSize',24)
    %legend('90°','180°','270°','360°')

    subplot(223)
    j = N_s/4;
    s_pos = 100*(j/N_s)*L;
    k0 = 1;
    k1 = N_th/4+1;
    k2 = N_th/2+1;
    k3 = 3*N_th/4+1;
    k4 = N_th+1;
    for i = 1:mov_CC_P
        P_k0(i) = CC_P(j,k0,i)*1.05;
        P_k1(i) = CC_P(j,k1,i)*1.05;
        P_k2(i) = CC_P(j,k2,i)*1.05;
        P_k3(i) = CC_P(j,k3,i)*1.05;
        P_k4(i) = CC_P(j,k4,i)*1.05;
    end
    plot(axis_T,P_k0,'b','Linewidth',4)
    hold on
    plot(axis_T,P_k1,'k','Linewidth',4)
    hold on
    plot(axis_T,P_k2,'r--','Linewidth',4)
    hold on
    plot(axis_T,P_k3,'g--','Linewidth',4)
    hold on
    %axis([0 2 -0.2 0.8])
    %hold on
    xlabel('Time (s)','FontSize',24)
    ylabel('Pressure (kPa)','FontSize',24)
	set(gca,'FontSize',24)
    %legend('90°','180°','270°','360°')

    sgtitle(['Profiles at s = ', num2str(s_pos),' cm'])
    hold on
%% %% Strain and Stress

	N_s = 200+4;
    N_th = 180;
    P = load('Example/P.dat');
    PT = load('Example/PT.dat');
    R_I = load('Example/R.dat');
    R_0 = load('Example/IC_Ro.dat');
    G_0 = load('Example/IC_Go.dat');
    for j = 1:size(R_I,1)
        R_E(j,1) = R_I(j,1) + (0.3)/(1000.0);
        R(j,1) = (R_I(j,1) + R_E(j,1))/2;
    end
    mov_R = size(R,1)/(N_s)/(N_th+1)

    R_0_Fixed = 1.0/100; %metros
    grosor = 0.03; % milimetros
    G_o = 40000;
    beta = 2.0;
    len = 2000;
    R_I_New = linspace(0.95/100.0,1.25/100.0,len);
    R_I_New = reshape(R_I_New,len,1);
    for i=1:len
        R_E_New(i,1) = R_I_New(i,1) + grosor/1000.0;
        R_prom(i,1) = 0.5*(R_I_New(i,1) + R_E_New(i,1));
        Strain_New(i,1) = R_prom(i,1)/R_0_Fixed;
        Stress_New(i,1) = 2.0*(G_o*((R_prom(i,1)/R_0_Fixed)^beta - 1) + 13.33)*(R_E_New(i,1)*R_I_New(i,1))^2/(R_E_New(i,1)^2 - R_I_New(i,1)^2)/R_prom(i,1)^2;
    end
    plot(Strain_New(:,1),Stress_New(:,1))
    hold on
    xlabel('Strain','FontSize',24)
    hold on
    ylabel('Stress','FontSize',24)
    hold on
    title('Strain vs. Stress, G_o = 1','FontSize',24)
    axis([1 1.2 0 2*10^7])

    for i = 1:mov_R
        for j = 1:N_s*(N_th+1)
            k = j+(i-1)*N_s*(N_th+1);
            %Strain(k,1) = R(k,1)/R_0_Fixed;
            Strain(k,1) = R(k,1)/R_0(j,1);
            Stress(k,1) = 2*(P(k,1) + 13.33)*(R_E(k,1)*R_I(k,1))^2/(R_E(k,1)^2 - R_I(k,1)^2)*(1/R(k,1)^2);
            %Stress(k,1) = 2*(G_o*((R/R_0)^beta-1));
        end
    end
 	R = reshape(R,N_s,N_th+1,mov_R);
 	P = reshape(P,N_s,N_th+1,mov_R);
	R_I = reshape(R_I,N_s,N_th+1,mov_R);
	R_E = reshape(R_E,N_s,N_th+1,mov_R);
    Strain = reshape(Strain,N_s,N_th+1,mov_R);
    Stress = reshape(Stress,N_s,N_th+1,mov_R);
    Data = load('Extracted.csv')
    Movie = VideoWriter('Example/Ex201/Videos/Stress-Strain_Horiz.avi');
	Movie.Quality = 100;
	Movie.FrameRate = 3;
    angle1 = 1;
    angle2 = 46;
    angle3 = 90;
    angle4 = 136;
    for mov=1:mov_R
        hold off
%        plot(Strain(3:1:N_s-2,angle,mov),Stress(3:1:N_s-2,angle,mov))
%        hold on
        plot(Strain(3:1:N_s-2,angle1,mov),Stress(3:1:N_s-2,angle1,mov),'.')
        hold on
        plot(Strain(3:1:N_s-2,angle2,mov),Stress(3:1:N_s-2,angle2,mov),'.')
        hold on
        plot(Strain(3:1:N_s-2,angle3,mov),Stress(3:1:N_s-2,angle3,mov),'.')
        hold on
        plot(Strain(3:1:N_s-2,angle4,mov),Stress(3:1:N_s-2,angle4,mov),'.')
        hold on
        plot(Data(:,1),Data(:,2)./10000.0,'.')
        hold on
        %axis([0.95 1.25 50 800])
        %hold on
        xlabel('Strain','FontSize',24)
        hold on
        ylabel('Stress','FontSize',24)
        hold on
        time = (mov-1)*Tmov;
        legend('0°','90°','180°','270°','DogData')
        title(['Strain vs. Stress at t = ', num2str(time),' seconds'])
        mov
	  	open(Movie)
	 	F = getframe(figure(1));
	 	writeVideo(Movie,F);
        pause(0.1);
    end
    close(Movie);


    title(['Stress-Strain at 1 second'])
    hold on
    angle = 25;
    for angle = 1:61
        subplot(231)
            subtitle(['Plot at 0 seconds'])
            hold on
            plot(Strain(:,angle,1),Stress(:,angle,1))
            hold on
            xlabel('Strain','FontSize',16)
            ylabel('Stress','FontSize',16)
        subplot(232)
            subtitle(['Plot at 0.2 seconds'])
            hold on
            plot(Strain(:,angle,21),Stress(:,angle,21))
            hold on
            xlabel('Strain','FontSize',16)
            ylabel('Stress','FontSize',16)
        subplot(233)
            subtitle(['Plot at 0.4 seconds'])
            hold on
            plot(Strain(:,angle,41),Stress(:,angle,41))
            hold on
            xlabel('Strain','FontSize',16)
            ylabel('Stress','FontSize',16)
        subplot(234)
            subtitle(['Plot at 0.6 seconds'])
            hold on
            plot(Strain(:,angle,61),Stress(:,angle,61))
            hold on
            xlabel('Strain','FontSize',16)
            ylabel('Stress','FontSize',16)
        subplot(235)
            subtitle(['Plot at 0.8 seconds'])
            hold on
            plot(Strain(:,angle,81),Stress(:,angle,81))
            hold on
            xlabel('Strain','FontSize',16)
            ylabel('Stress','FontSize',16)
        subplot(236)
            subtitle(['Plot at 1 seconds'])
            hold on
            plot(Strain(:,angle,101),Stress(:,angle,101))
            hold on
            xlabel('Strain','FontSize',16)
            ylabel('Stress','FontSize',16)
        hold off
        angle
        pause(0.1)
    end
    print('Example/Ex201/Videos/Strain-Stress','-depsc',figure(1))

    
%% LOAD PARAMETERS

    N_s = 100+4;
    N_th = 60;
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
    
    
    PT = load('Example/Ex96/Files/PT.dat');
    P = load('Example/Ex96/Files/P.dat');
    R_I = load('Example/Ex96/Files/R.dat');
    R_0 = load('Example/Ex96/IC/IC_Ro.dat');
    G_0 = load('Example/Ex96/IC/IC_Go.dat');
    for j = 1:size(R_I,1)
        R_E(j,1) = R_I(j,1) + (0.03)/(1000.0);
        R(j,1) = (R_I(j,1) + R_E(j,1))/2;
    end
    mov_R = size(R,1)/(N_s)/(N_th+1)
    
    
    for i = 1:mov_R
        for j = 1:N_s*(N_th+1)
            k = j+(i-1)*N_s*(N_th+1);
            Strain_2D(k,1) = R(k,1)/R_0(j,1);
            Stress_2D(k,1) = 2*(1.05*P(k,1) + P0)*(R_E(k,1)*R_I(k,1))^2/(R_E(k,1)^2 - R_I(k,1)^2)*(1/R(k,1)^2);
        end
    end
 	R = reshape(R,N_s,N_th+1,mov_R);
 	P = reshape(P,N_s,N_th+1,mov_R);
	R_I = reshape(R_I,N_s,N_th+1,mov_R);
	R_E = reshape(R_E,N_s,N_th+1,mov_R);
    G_0 = reshape(G_0,N_s,N_th+1,2);
    Strain_2D = reshape(Strain_2D,N_s,N_th+1,mov_R);
    Stress_2D = reshape(Stress_2D,N_s,N_th+1,mov_R);
    Data = load('Extracted.csv')

    Movie = VideoWriter('Example/Ex96/Videos/Stress-Strain_Horiz.avi');
	Movie.Quality = 100;
	Movie.FrameRate = 3;
    angle1 = 1;
    angle2 = 16;
    angle3 = 31;
    angle4 = 46;
    for mov=1:mov_R
        hold off
%        plot(Strain(3:1:N_s-2,angle,mov),Stress(3:1:N_s-2,angle,mov))
%        hold on
        plot(Strain_2D(3:1:N_s-2,angle1,mov),Stress_2D(3:1:N_s-2,angle1,mov),'.')
        hold on
        %plot(Strain_2D(3:1:N_s-2,angle2,mov),Stress_2D(3:1:N_s-2,angle2,mov),'.')
        %hold on
        %plot(Strain_2D(3:1:N_s-2,angle3,mov),Stress_2D(3:1:N_s-2,angle3,mov),'.')
        %hold on
        %plot(Strain_2D(3:1:N_s-2,angle4,mov),Stress_2D(3:1:N_s-2,angle4,mov),'.')
        %hold on
        plot(Data(:,1),Data(:,2)./10000.0,'.')
        hold on
        %axis([0.95 1.25 50 1500])
        %hold on
        xlabel('Strain','FontSize',24)
        hold on
        ylabel('Stress','FontSize',24)
        hold on
        time = (mov-1)*Tmov;
        legend('0°','90°','180°','270°')
        title(['Strain vs. Stress at t = ', num2str(time),' seconds'])
        mov
	  	open(Movie)
	 	F = getframe(figure(1));
	 	writeVideo(Movie,F);
        pause(0.1);
    end
    close(Movie);
    
    
    P0 = 1;
    for i=1:2000
    Stress_P0(1,i) = 2*P0/(RRI(1,i).^2 - RRE(1,i).^2);
    end
    
    
    
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

%% valores que aproximan los datos
    clf
    figure(1)
    Data = load('Extracted.csv')
    Strain = linspace(0.8,1.5,2000);
    h = 0.03; % milimetros
    R0 = 1.0; % centimetros
    beta = 0.25; % adimensional
    G0 = 5.0; % kilopascales
    P0 = 0.17; %kilopascales
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
    plot(Data(:,1),Data(:,2)./10000.0,'.')
    hold on
    axis([1.09 1.23 50 200])
    hold on
    print('Strain-Stress_Approach','-depsc',figure(1))

    
%% contraste con 2D modelo    
end
