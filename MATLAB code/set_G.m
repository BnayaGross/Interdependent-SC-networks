function [G, state] = set_G(W,L,N,R_ij,Rsmall,Ic_T)
G = zeros(N);
state = zeros(N);
%supernodes

G(1,1) = 0;
G(N,N) = 0;
for i = 1:L
    % off diagonal

    Vij = abs(W(i+1) - W(1));
    if Vij > R_ij(1,i+1)*Ic_T(1,i+1)
        G(1,i+1) = -1/R_ij(1,i+1);
        state(1,i+1) = 1;
    elseif Vij < Rsmall*Ic_T(1,i+1)
        G(1,i+1) = -1/Rsmall;
        state(1,i+1) = 3;
    else
        G(1,i+1) = -Ic_T(1,i+1)/Vij;
        state(1,i+1) = 2;
    end

    G(i+1,1) = G(1,i+1);
    state(i+1,1) = state(1,i+1);

    %diagonal
    G(1,1) = G(1,1) - G(1,i+1);


    Vij = abs(W(N) - W(N - i));
    if Vij > R_ij(N,N - i)*Ic_T(N,N - i)
        G(N - i,N) = -1/R_ij(N,N - i);
        state(N-i, N) = 1;
    elseif Vij < Rsmall*Ic_T(N,N - i)
        G(N - i,N) = -1/Rsmall;
        state(N-i, N) = 3;
    else
        G(N - i,N) = -Ic_T(N,N - i)/Vij;
        state(N-i, N) = 2;
    end

    G(N,N - i) = G(N - i, N);
    state(N,N - i) = state(N - i, N);  

    %diagonal
    G(N,N) = G(N,N) - G(N,N-i);
end

for i = 2:N-1

    %off diagonal

    % down
    if mod(i-1,L) ~= 0
        Vij = abs(W(i+1) - W(i));
        if Vij > R_ij(i,i+1)*Ic_T(i,i+1)
            G(i,i+1) = -1/R_ij(i,i+1);
            state(i,i+1) = 1;
        elseif Vij < Rsmall*Ic_T(i,i+1)
            G(i,i+1) = -1/Rsmall;
            state(i,i+1) = 3;
        else
            G(i,i+1) = -Ic_T(i,i+1)/Vij;
            state(i,i+1) = 2;
        end

        G(i+1,i) = G(i,i+1);
        state(i+1,i) = state(i,i+1);
    end

    %right
    if i < N - L
        Vij = abs(W(i+L) - W(i));
        if Vij > R_ij(i,i+L)*Ic_T(i,i+L)
            G(i,i+L) = -1/R_ij(i,i+L);
            state(i,i+L) = 1;
        elseif Vij < Rsmall*Ic_T(i,i+L)
            G(i,i+L) = -1/Rsmall;
            state(i,i+L) = 3;
        else
            G(i,i+L) = -Ic_T(i,i+L)/Vij;
            state(i,i+L) = 2;
        end

        G(i+L,i) = G(i,i+L);  
        state(i+L,i) = state(i,i+L);  
    end


    %diagonal

    G(i,i) = 0;


    %down
    if mod(i-1,L) ~= 0
        G(i,i) = G(i,i) - G(i,i+1);
    end

    %up
    if mod(i-2,L) ~= 0
        G(i,i) = G(i,i) - G(i,i-1);
    end

    %left
    if i > L + 1
        G(i,i) = G(i,i) - G(i,i-L);
    else
        G(i,i) = G(i,i) - G(i,1);
    end

    %right
    if i < N-L
        G(i,i) = G(i,i) - G(i,i+L);
    else
        G(i,i) = G(i,i) - G(i,N);
    end
end

