%% Capsula ARD - Esercizio rientro orbitale
% Calcolo la traiettoria di rientro e i profili di decelerazione e flusso
% termico della capsula ARD durante il rientro orbitale in atmosfera terrestre
% Codice accoppiato alla funzione missionreturn.m
clc
clear
close all

% Legge il file che contiene i dati dell'atmosfera standard
A = dlmread('stdatm.dat');  

% Dati Importati della stdatm
Hd = A(:,1);     % Altezza al variare della quota
R0 = 6.37*10^6;  % Raggio Terrrestre
Rd = R0 + Hd;    % Raggio dell'orbita 
Td  = A(:,2);    % Temperatura al variare della quota
rhod = A(:,4);   % Densita al variare della quota

% Parametri da passare alla funzione integranda missionreturn.m
Mass = 2800;
Diam = 2.8;
Area_f = pi*(Diam/2)^2;
C_D = 1;
Balistic_coeff = Mass/(Area_f*C_D);
Bb   = [100,Balistic_coeff,700]; % Coefficiente Balistico [kg/m2]
L_D = 0;                         % Efficienza aerodinamica (L/D) (No Lift)
Rc  = 1;                         % Curvature radius [m]
eps = 0.8;                       % Surface emissivity

% CONDIZIONI INIZIALI 
% Ipotizzati valori costandi di Stanton e Coeff. Balistico (CD = 1, St = 0.01)
H_i = 120;                       % Quota iniziale [km]
V_i = 7500;                      % Velocita iniziale [m/s]
gamma_i = [1.5,5,9];             % Angolo di rientro iniziale / flight path angle [deg]
St = 0.01;                       % Stanton 
R=287;
Cp = 1.4*R/(0.4);

% Analisi due casi
caso = 'gamma'; % or 'B'
switch caso
    case 'gamma'  % Fisso il coeff. balistico e vario gamma iniziale
 % Fisso il coeff. balistico e vario gamma
        B = Bb(2) % B = 454
        Br = 1./B;
        for i = 1:3 % gamma_i = [1.5,5,9]; 
[V1, gamma1, H1, time1, rho1, T1, g1]=missionreturn(B, L_D, Rc, eps, H_i, V_i, gamma_i(i));
% Risultati forniti dalla funzione missionreturn.m
% La funzione calcola l'andamento della velocita e dell'angolo gamma
% mediante il metodo di integrazione di Eulero
            l = length(V1);
            V(1:l,i) = V1;
            H(1:l,i) = H1;
            gamma(1:l,i) = gamma1; % Assegnato il gamma iniziale, durante la missione esso variera
            time(1:l,i) = time1;
            rho(1:l,i) = rho1;
            T(1:l,i) = T1;
            g(1:l,i) = g1;  

            % Calcolo Profili di decelerazione, flussi termici ecc.
            acc(1:l,i)  = (-0.5*rho1.*V1.^2.*Br+g1.*sin(gamma1)); % Decelerazione
            M(1:l,i)  = V1./sqrt(1.4*287*T1);  % Mach
            p0(1:l,i) = rho1.*V1.^2; % pressione di ristagno
            T_adiabatic_wall = 900;
            q(1:l,i) =  rho1.*V1*St.*((V1.^2/2) + Cp*(T1 - T_adiabatic_wall));    % Flusso di calore
        end

        %% ANALITICA 
        % Per confrontare i risultati numerici con quelli analitici
        % In questo caso si suppone che:
        % Gamma è costante per tutta la missione
        % Non vi è la forza peso g=0
        % Modello di atmosfera isoterma
        
        gamma_in_rad = convang(gamma_i,'deg','rad'); % Studio il problema al variare del gamma iniziale
        N = 570;
        H_an = linspace(0,H_i*1000,N); % Vettore delle altezze
        H_an = fliplr(H_an);
        % Modello di atmosfera isoterma
        beta = 1/8500; 
        z = beta*H_an;
        rho_0 = 1.225;  
        
        for i=1:3 % % gamma_i = [1.5,5,9]; 
            
            alpha = (1\2)*(rho_0/(beta*B*sin(gamma_in_rad(i)))); % Parametro adimensionale
        
 % Calcolo incognite analiticamente al variare della quota e del gamma iniziale
            for k = 1:N
                V_an(k,i) = V_i*exp(-alpha*exp(-z(k)));
                acc_an(k,i) = -alpha*beta*sin(gamma_in_rad(i))*(V_an(k,i)^(2))*exp(-z(k));
                q_punto_an(k,i) = rho_0*((V_an(k,i)^3)/2)*St*exp(-z(k));
                a_suono = 305;
                M_an(k,i) = V_an(k,i)/a_suono;
            end

        end
        
%% plot
        figure
        plot(V(:,1),H(:,1)/1000,'k','linewidth',1.5);
        hold on
        plot(V(:,2),H(:,2)/1000,'-.k','linewidth',1.5);
        hold on;
        plot(V(:,3),H(:,3)/1000,':k','linewidth',1.5);
        hold on;
        plot(V_an(:,1),H_an(:)/1000,'r','linewidth',1.5);
        hold on;
        plot(V_an(:,2),H_an(:)/1000,'-.r','linewidth',1.5);
        hold on;
        plot(V_an(:,3),H_an(:)/1000,':r','linewidth',1.5);
        grid on
        xlabel('Velocities [m/s]')
        ylabel('H [Km]')
        set(gcf,'color','w')
        legend('\gamma_i = 1.5','\gamma_i = 5', '\gamma_i = 9')

        figure
        plot(M(:,1),H(:,1)/1000,'k','linewidth',1.5);
        hold on
        plot(M(:,2),H(:,2)/1000,'-.k','linewidth',1.5);
        hold on;
        plot(M(:,3),H(:,3)/1000,':k','linewidth',1.5);
        hold on;
        plot(M_an(:,1),H_an(:)/1000,'r','linewidth',1.5);
        hold on;
        plot(M_an(:,2),H_an(:)/1000,'-.r','linewidth',1.5);
        hold on;
        plot(M_an(:,3),H_an(:)/1000,':r','linewidth',1.5);
        grid on
        xlabel('Mach')
        ylabel('H [Km]')
        set(gcf,'color','w')
        legend('\gamma_i = 1.5','\gamma_i = 5', '\gamma_i = 9')        
        

        figure
        idx2 = time(:,2)>0;
        idx3 = time(:,3)>0;
        plot(time(:,1),H(:,1)/1000,'k','linewidth',1.5);
        hold on
        plot(time(idx2,2),H(idx2,2)/1000,'-.k','linewidth',1.5)
        plot(time(idx3,3),H(idx3,3)/1000,':k','linewidth',1.5)
        grid on
        xlabel('Time [s]')
        ylabel('H [Km]')
        legend('\gamma_i = 1.5','\gamma_i = 5', '\gamma_i = 9')

        figure
        idx = time>0;
        plot(time(idx(:,1),1),V(idx(:,1),1),'--k','linewidth',1.5)
        hold on
        plot(time(idx(:,2),2),V(idx(:,2),2),'-.k','linewidth',1.5)
        plot(time(idx(:,3),3),V(idx(:,3),3),':k','linewidth',1.5)
        xlabel('Time [s]')
        ylabel('Velocities [m/s]')
        set(gcf,'color','w')
        legend('\gamma_i = 1.5','\gamma_i = 5', '\gamma_i = 9')

        figure
        plot(-gamma(:,1)*180/pi,H(:,1)/1000,'k','linewidth',1.5)
        hold on
        plot(-gamma(:,2)*180/pi,H(:,2)/1000,'-.k','linewidth',1.5)
        plot(-gamma(:,3)*180/pi,H(:,3)/1000,':k','linewidth',1.5)
        hold on
        grid on
        xlabel('\gamma [° ]')
        ylabel('H [Km]')
        set(gcf,'color','w')
        legend('\gamma_i = 1.5','\gamma_i = 5', '\gamma_i = 9')

        figure
        plot(acc(:,1)/9.81,H(:,1)/1000,'k','linewidth',1.5)
        hold on
        plot(acc(:,2)/9.81,H(:,2)/1000,'-.k','linewidth',1.5)
        hold on
        plot(acc(:,3)/9.81,H(:,3)/1000,':k','linewidth',1.5)
        hold on
        plot(acc_an(:,1),H_an(:)/1000,'r','linewidth',1.5);
        hold on;
        plot(acc_an(:,2),H_an(:)/1000,'-.r','linewidth',1.5);
        hold on;
        plot(acc_an(:,3),H_an(:)/1000,':r','linewidth',1.5);
        grid on
        xlabel('a/go')
        ylabel('H [Km]')
        set(gcf,'color','w')
        legend('\gamma_i = 1.5','\gamma_i = 5', '\gamma_i = 9')

        figure
        plot(q(:,1),H(:,1)/1000,'k','linewidth',1.5)
        hold on
        plot(q(:,2),H(:,2)/1000,'-.k','linewidth',1.5)
        hold on
        plot(q(:,3),H(:,3)/1000,':k','linewidth',1.5)
        hold on
        plot(q_punto_an(:,1),H_an(:)/1000,'r','linewidth',1.5);
        hold on;
        plot(q_punto_an(:,2),H_an(:)/1000,'-.r','linewidth',1.5);
        hold on;
        plot(q_punto_an(:,3),H_an(:)/1000,':r','linewidth',1.5);        
        grid on
        xlabel('q [W/m2]')
        ylabel('H [Km]')
        set(gcf,'color','w')
        legend('\gamma_i = 1.5','\gamma_i = 5', '\gamma_i = 9')

 case 'B' % Fisso gamma = 1.5 e vario il coeff Balistico 
        gamma_i = gamma_i(1)

        for i = 1:3
            B = Bb;
            Br = 1./B(i);
            [V1, gamma1, H1, time1, rho1, T1, g1]=missionreturn(B(i), L_D, Rc, eps, H_i, V_i, gamma_i);

            l = length(V1);
            V(1:l,i) = V1;
            H(1:l,i) = H1;
            gamma(1:l,i) = gamma1;
            time(1:l,i) = time1;
            rho(1:l,i) = rho1;
            T(1:l,i) = T1;
            g(1:l,i) = g1;

            acc(1:l,i)  = (-0.5*rho1.*V1.^2.*Br+g1.*sin(gamma1))/9.81;
            M(1:l,i)  = V1./sqrt(1.4*287*T1);                % Mach
            %q  = (1.83e-4)*sqrt(rho/Rc).*V.^3;      % flusso di calore con formula di Tauber
            p0(1:l,i) = rho1.*V1.^2;                         % pressione di ristagno
            q(1:l,i) =  rho1.*V1*St.*(V1.^2/2);% +Cp*(T - 900));
        end

        figure
        plot(V(:,1),H(:,1)/1000,'k','linewidth',1.5);
        hold on
        plot(V(:,2),H(:,2)/1000,'-.k','linewidth',1.5);
        plot(V(:,3),H(:,3)/1000,':k','linewidth',1.5);
        grid on
        xlabel('Velocities [m/s]')
        ylabel('H [Km]')
        set(gcf,'color','w')
        legend('B = 100','B = 450', 'B = 700')

        figure
        idx2 = time(:,2)>0;
        idx3 = time(:,3)>0;
        plot(time(:,1),H(:,1)/1000,'k','linewidth',1.5);
        hold on
        plot(time(idx2,2),H(idx2,2)/1000,'-.k','linewidth',1.5)
        plot(time(idx3,3),H(idx3,3)/1000,':k','linewidth',1.5)
        grid on
        xlabel('Time [s]')
        ylabel('H [Km]')
        legend('B = 100','B = 450', 'B = 700')


        figure
        idx = time>0;
        plot(time(idx(:,1),1),V(idx(:,1),1),'--k','linewidth',1.5)
        hold on
        plot(time(idx(:,2),2),V(idx(:,2),2),'-.k','linewidth',1.5)
        plot(time(idx(:,3),3),V(idx(:,3),3),':k','linewidth',1.5)
        xlabel('Time [s]')
        ylabel('Velocities [m/s]')
        set(gcf,'color','w')
        legend('B = 100','B = 450', 'B = 700')


        figure
        plot(gamma(:,1)*180/pi,H(:,1)/1000,'k','linewidth',1.5)
        hold on
        plot(gamma(:,2)*180/pi,H(:,2)/1000,'-.k','linewidth',1.5)
        plot(gamma(:,3)*180/pi,H(:,3)/1000,':k','linewidth',1.5)
        hold on
        grid on
        xlabel('\gamma [° ]')
        ylabel('H [Km]')
        set(gcf,'color','w')
        legend('B = 100','B = 450', 'B = 700')


        figure
        plot(acc(:,1),H(:,1)/1000,'k','linewidth',1.5)
        hold on
        plot(acc(:,2),H(:,2)/1000,'-.k','linewidth',1.5)
        plot(acc(:,3),H(:,3)/1000,':k','linewidth',1.5)
        hold on
        grid on
        xlabel('a [m/s2]')
        ylabel('H [Km]')
        set(gcf,'color','w')
        legend('B = 100','B = 450', 'B = 700')


        figure
        plot(q(:,1),H(:,1),'k','linewidth',1.5)
        hold on
        plot(q(:,2),H(:,2),'-.k','linewidth',1.5)
        plot(q(:,3),H(:,3),':k','linewidth',1.5)
        hold on
        grid on
        xlabel('q [W/m2]')
        ylabel('H [Km]')
        set(gcf,'color','w')
        legend('B = 100','B = 450', 'B = 700')
end
