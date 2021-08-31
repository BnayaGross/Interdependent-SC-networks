%% Interdependent SC networks

% Here we will simulate the case of two coupled layers and then compare it
% to the isolated case on the exact same structure.



tests = 1:10;
for idx_t = tests

    %% parameters
    L = 31;
    N = L*L + 2;
    R0 = 500;
    Rsmall = 10^-5;
    Ico1 = 58*10^-6;
    Ico2 = 58*10^-6;
    dT = 10^-3;
    Ti = 0.001;
    Tf = 20;
    sigma = 0.1;
    Ib = 24*10^-6;

    a1 = 5*10^5; %%heating inside layer 1
    a2 = 5*10^5; %%heating inside layer 2
    a12 = 3*10^6; %%heating from layer 2 to layer 1
    a21 = 3*10^6; %%heating from layer 1 to layer 2

    C = 1;
    gamma = 2;
    movie_flag = 0;


    %% initialization
    R1 = 0;
    R2 = 0;

    W1 = zeros(N,1);
    W1(N) = 0;
    W2 = zeros(N,1);
    W2(N) = 0;

    Ib_vec1 = zeros(N-1,1);
    Ib_vec1(1) = Ib;
    Ib_vec2 = zeros(N-1,1);
    Ib_vec2(1) = Ib;

    R_heating_coupled_vec1 = [];
    R_heating_coupled_vec2 = [];
    T_heating_coupled_vec = [];
    iter_heating_coupled_vec = [];

    [Ic_zero1, Aij1] = set_Ic_zero(L,N, sigma);  %% comment this if you want to load
    [Ic_zero2, Aij2] = set_Ic_zero(L,N, sigma);  %% comment this if you want to load

    
    %Ic_zero1 = Ic_zero2; %% uncomment if you want identical layers
    %Aij1 = Aij2;
    
    %save setup3_two_layers_sigma0.1.mat Ic_zero1 Aij1 Ic_zero2 Aij2;    %% uncomment what you need
    %load setup_two_layers_sigma0.1.mat

    %R_ij1 = C./(Ico*Ic_zero1);
    %R_ij2 = C./(Ico*Ic_zero2);
    R_ij1 = R0*Aij1;    %%if you want to go back to fixed value of R0 uncomment this
    R_ij2 = R0*Aij2;    %%if you want to go back to fixed value of R0 uncomment this

    
    Tc1 = 100*Ico1*Ic_zero1*R0;
    Tc2 = 100*Ico2*Ic_zero2*R0;

    if movie_flag
        v = VideoWriter('newfile.avi');
        open(v);
    end

    flag1 = 0;
    flag2 = 0;

    %% heating coupled
    for T = Ti:dT:Tf

        epsilon_vec = [];

        for iter = 1:1000000


            Teff1 = T + (a1*R1 + a12*R2)*Ib^2;
            Teff2 = T + (a21*R1 + a2*R2)*Ib^2;


            if flag1 == 0
                T_temp1 = (1 - Teff1./Tc1);
                Ic_T1 = Ico1*Ic_zero1.*(heaviside(T_temp1).*T_temp1).^gamma;
                [G1, state1] = set_G(W1,L,N,R_ij1,Rsmall,Ic_T1);
                G21 = G1(1:N-1, 1:N-1);
                W_next1 = zeros(N,1);
                W_next1(1:N-1) = G21\Ib_vec1;
                W_next1(N) = 0;
                R1 = W_next1(1)/Ib;
            end


            if flag2 == 0
                T_temp2 = (1 - Teff2./Tc2);
                Ic_T2 = Ico2*Ic_zero2.*(heaviside(T_temp2).*T_temp2).^gamma;
                [G2, state2] = set_G(W2,L,N,R_ij2,Rsmall,Ic_T2);
                G22 = G2(1:N-1, 1:N-1);
                W_next2 = zeros(N,1);
                W_next2(1:N-1) = G22\Ib_vec2;
                W_next2(N) = 0;
                R2 = W_next2(1)/Ib;
            end


            if R1 > 1.03*R0 && flag1 == 0
                flag1 = 1;
                Tflag1 = T;
            end

            if R2 > 1.03*R0
                flag2 = 1;
            end

            if flag1 && flag2
                break
            end



            if any(isnan(W_next1)) || any(isnan(W_next2))
                break;
            end


            epsilon = (norm(W_next1 - W1))/(norm(W_next1)) + (norm(W_next2 - W2))/(norm(W_next2)); 
            epsilon_vec = [epsilon_vec epsilon];



            [T, epsilon, iter]
            R1

            W_pre1 = W1;
            W1 = W_next1;
            W_pre2 = W2;
            W2 = W_next2;

            if epsilon < 10^-5
                break;
            end
        end

        T_heating_coupled_vec = [T_heating_coupled_vec T];
        R_heating_coupled_vec1 = [R_heating_coupled_vec1 R1];
        R_heating_coupled_vec2 = [R_heating_coupled_vec2 R2];
        iter_heating_coupled_vec = [iter_heating_coupled_vec iter];


        %% add frame to the movie
        if movie_flag
            [frame] = set_frame_for_movie_two_layers(L,N,Ti,4.2,R0, state1,  state2, T_heating_vec, R_heating_vec1, R_heating_vec2, [],[],[]);
            writeVideo(v , frame);
            close;
        end


        if flag1 && flag2
            break
        end
    end

    %% cooling coupled
    R1 = R_heating_coupled_vec1(end);
    R2 = R_heating_coupled_vec2(end);
    new_Ti = T;
    R_cooling_coupled_vec1 = [];
    R_cooling_coupled_vec2 = [];
    T_cooling_coupled_vec = [];
    iter_cooling_coupled_vec = [];

    for T = new_Ti:-dT:Ti
        T
        epsilon_vec = [];

        for iter = 1:1000000

            if T <= Tflag1
                Teff1 = T + (a1*R1 + a12*R2)*Ib^2;
                T_temp1 = (1 - Teff1./Tc1);
                Ic_T1 = Ico1*Ic_zero1.*(heaviside(T_temp1).*T_temp1).^gamma; 
                [G1, state1] = set_G(W1,L,N,R_ij1,Rsmall,Ic_T1);
                G21 = G1(1:N-1, 1:N-1);
                W_next1 = zeros(N,1);
                W_next1(1:N-1) = G21\Ib_vec1;
                W_next1(N) = 0;
                R1 = W_next1(1)/Ib;
            end

            Teff2 = T + (a21*R1 + a2*R2)*Ib^2;
            T_temp2 = (1 - Teff2./Tc2);
            Ic_T2 = Ico2*Ic_zero2.*(heaviside(T_temp2).*T_temp2).^gamma; 
            [G2, state2] = set_G(W2,L,N,R_ij2,Rsmall,Ic_T2);
            G22 = G2(1:N-1, 1:N-1);
            W_next2 = zeros(N,1);
            W_next2(1:N-1) = G22\Ib_vec2;
            W_next2(N) = 0;
            R2 = W_next2(1)/Ib;

            epsilon = (norm(W_next1 - W1))/(norm(W_next1)) + (norm(W_next2 - W2))/(norm(W_next2)); 
            epsilon_vec = [epsilon_vec epsilon];

            [T epsilon]
            R1



            W_pre1 = W1;
            W1 = W_next1;
            W_pre2 = W2;
            W2 = W_next2;

            if epsilon < 10^-5
                break;
            end
        end

        T_cooling_coupled_vec = [T_cooling_coupled_vec T];
        R_cooling_coupled_vec1 = [R_cooling_coupled_vec1 R1];
        R_cooling_coupled_vec2 = [R_cooling_coupled_vec2 R2];
        iter_cooling_coupled_vec = [iter_cooling_coupled_vec iter];


        %% add frame to the movie
        if movie_flag
            [frame] = set_frame_for_movie_two_layers(L,N,Ti,4.2,R0, state1,  state2, T_heating_vec, R_heating_vec1, R_heating_vec2, T_cooling_vec, R_cooling_vec1, R_cooling_vec2);
            writeVideo(v , frame);
            close;
        end


        if R1 < 10^-3 && R2 < 10^-3
            break
        end
    end

    if movie_flag
        close(v);
    end

    figure(); hold on;
    % plot(T_heating_coupled_vec, R_heating_coupled_vec1/max(R_heating_coupled_vec1), '-ro');
    % plot(T_heating_coupled_vec, R_heating_coupled_vec2/max(R_heating_coupled_vec2), '-bo');
    % plot(T_cooling_coupled_vec, R_cooling_coupled_vec1/max(R_cooling_coupled_vec1), '-rs');
    % plot(T_cooling_coupled_vec, R_cooling_coupled_vec2/max(R_cooling_coupled_vec2), '-bs');
    plot(T_heating_coupled_vec, R_heating_coupled_vec1, '-ro');
    plot(T_heating_coupled_vec, R_heating_coupled_vec2, '-bo');
    plot(T_cooling_coupled_vec, R_cooling_coupled_vec1, '-rs');
    plot(T_cooling_coupled_vec, R_cooling_coupled_vec2, '-bs');
    %xlim([4,4.6]);
    %ylim([0,1.1]);
    legend('heating layer1','heating layer2','cooling layer1','cooling layer2', 'Location', 'best');
    %title('\alpha = 1000, I_b = 1 \mu A, I_{c0} = 0.5 \mu A,R_0 = 1.5 \cdot 10^3');
    ylabel('R');
    xlabel('T(K)');


    R1 = 0;
    R2 = 0;

    W1 = zeros(N,1);
    W1(N) = 0;
    W2 = zeros(N,1);
    W2(N) = 0;

    R_heating_sep_vec1 = [];
    R_heating_sep_vec2 = [];
    T_heating_sep_vec = [];
    iter_heating_sep_vec = [];

    flag1 = 0;
    flag2 = 0;

    %% heating separated layers
    for T = Ti:dT:Tf
        T
        epsilon_vec = [];

        for iter = 1:100000

            Teff1 = T + a1*R1*Ib^2;
            Teff2 = T + a2*R2*Ib^2;

            if flag1 == 0
                T_temp1 = (1 - Teff1./Tc1);
                Ic_T1 = Ico1*Ic_zero1.*(heaviside(T_temp1).*T_temp1).^gamma;
                [G1, state1] = set_G(W1,L,N,R_ij1,Rsmall,Ic_T1);
                G21 = G1(1:N-1, 1:N-1);
                W_next1 = zeros(N,1);
                W_next1(1:N-1) = G21\Ib_vec1;
                W_next1(N) = 0;
                epsilon1 = (norm(W_next1 - W1))/(norm(W_next1));
                R1 = W_next1(1)/Ib;
            end

            if flag2 == 0
                T_temp2 = (1 - Teff2./Tc2);
                Ic_T2 = Ico2*Ic_zero2.*(heaviside(T_temp2).*T_temp2).^gamma; 
                [G2, state2] = set_G(W2,L,N,R_ij2,Rsmall,Ic_T2);
                G22 = G2(1:N-1, 1:N-1);
                W_next2 = zeros(N,1);
                W_next2(1:N-1) = G22\Ib_vec2;
                W_next2(N) = 0;
                epsilon2 = (norm(W_next2 - W2))/(norm(W_next2));
                R2 = W_next2(1)/Ib;
            end


            if any(isnan(W_next1)) || any(isnan(W_next2))
                break;
            end



            [T, epsilon1,epsilon2 , iter]
            R1

            if R1 > 1.03*R0 && flag1 == 0
                flag1 = 1;
                Tflag1 = T;
            end

            if R2 > 1.03*R0
                flag2 = 1;
            end

            if flag1 && flag2
                break
            end


            W_pre1 = W1;
            W1 = W_next1;
            W_pre2 = W2;
            W2 = W_next2;

            if epsilon1 < 10^-5 && epsilon2 < 10^-5
                break;
            end

            if flag1 && epsilon2 < 10^-5
                break;
            end

            if flag2 && epsilon1 < 10^-5
                break;
            end
        end

        T_heating_sep_vec = [T_heating_sep_vec T];
        R_heating_sep_vec1 = [R_heating_sep_vec1 R1];
        R_heating_sep_vec2 = [R_heating_sep_vec2 R2];
        iter_heating_sep_vec = [iter_heating_sep_vec iter];


        if flag1 && flag2
            break
        end
    end


    %% cooling separated layers
    R1 = R_heating_sep_vec1(end);
    R2 = R_heating_sep_vec2(end);
    new_Ti = T;
    R_cooling_sep_vec1 = [];
    R_cooling_sep_vec2 = [];
    T_cooling_sep_vec = [];
    iter_cooling_sep_vec = [];

    for T = new_Ti:-dT:Ti
        T
        epsilon_vec = [];

        for iter = 1:100000

            if T <= Tflag1
                Teff1 = T + a1*R1*Ib^2;
                T_temp1 = (1 - Teff1./Tc1);
                Ic_T1 = Ico1*Ic_zero1.*(heaviside(T_temp1).*T_temp1).^gamma; 
                [G1, state1] = set_G(W1,L,N,R_ij1,Rsmall,Ic_T1);
                G21 = G1(1:N-1, 1:N-1);
                W_next1 = zeros(N,1);
                W_next1(1:N-1) = G21\Ib_vec1;
                W_next1(N) = 0;
                R1 = W_next1(1)/Ib;
            end

            Teff2 = T + a2*R2*Ib^2;
            T_temp2 = (1 - Teff2./Tc2);
            Ic_T2 = Ico2*Ic_zero2.*(heaviside(T_temp2).*T_temp2).^gamma; 
            [G2, state2] = set_G(W2,L,N,R_ij2,Rsmall,Ic_T2);
            G22 = G2(1:N-1, 1:N-1);
            W_next2 = zeros(N,1);
            W_next2(1:N-1) = G22\Ib_vec2;
            W_next2(N) = 0;
            R2 = W_next2(1)/Ib;


            epsilon1 = (norm(W_next1 - W1))/(norm(W_next1));
            epsilon2 = (norm(W_next2 - W2))/(norm(W_next2));
           % epsilon_vec = [epsilon_vec epsilon];


            W_pre1 = W1;
            W1 = W_next1;
            W_pre2 = W2;
            W2 = W_next2;

            [T, epsilon1,epsilon2 , iter]
            R1

            if epsilon1 < 10^-5 && epsilon2 < 10^-5
                break;
            end
        end

        T_cooling_sep_vec = [T_cooling_sep_vec T];
        R_cooling_sep_vec1 = [R_cooling_sep_vec1 R1];
        R_cooling_sep_vec2 = [R_cooling_sep_vec2 R2];
        iter_cooling_sep_vec = [iter_cooling_sep_vec iter];



        if R1 < 10^-3 && R2 < 10^-3
            break
        end
    end


    figure(); hold on;
    % plot(T_heating_sep_vec, R_heating_sep_vec1/max(R_heating_sep_vec1), '-ro');
    % plot(T_heating_sep_vec, R_heating_sep_vec2/max(R_heating_sep_vec2), '-bo');
    % plot(T_cooling_sep_vec, R_cooling_sep_vec1/max(R_cooling_sep_vec1), '-rs');
    % plot(T_cooling_sep_vec, R_cooling_sep_vec2/max(R_cooling_sep_vec2), '-bs');
    plot(T_heating_sep_vec, R_heating_sep_vec1, '-ro');
    plot(T_heating_sep_vec, R_heating_sep_vec2, '-bo');
    plot(T_cooling_sep_vec, R_cooling_sep_vec1, '-rs');
    plot(T_cooling_sep_vec, R_cooling_sep_vec2, '-bs');

    %xlim([4,4.6]);
    %ylim([0,1.1]);
    legend('heating layer1','heating layer2','cooling layer1','cooling layer2', 'Location', 'best');

    %title('a11 = 2000, a12 = 130000, I_{co1} = 70 \muA, I_{co2} = 7 \muA, \sigma = 1, I_b = 24 \mu A, T_c = 2.8');
    ylabel('R_N');
    xlabel('T(k)');


    %filename = ['Ic148Ic258sigma0.1alpha3Ib' ,num2str(Ib*(10^6)),'_h_test', num2str(idx_t), '.mat'];
        filename = ['Ic148Ic258sigma0.1alpha3Ib' ,num2str(Ib*(10^6)),'_identical_test', num2str(idx_t), '.mat'];

    save(filename, 'T_heating_sep_vec','T_cooling_sep_vec','R_heating_sep_vec1','R_heating_sep_vec2','R_cooling_sep_vec1','R_cooling_sep_vec2','T_heating_coupled_vec','T_cooling_coupled_vec','R_heating_coupled_vec1','R_heating_coupled_vec2','R_cooling_coupled_vec1','R_cooling_coupled_vec2');
    clc;
    clear;
end

