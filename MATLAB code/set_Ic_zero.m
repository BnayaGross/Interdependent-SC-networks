function [Ic_zero, Aij] = set_Ic_zero(L,N, sigma)
Ic_zero = zeros(N);
Aij = zeros(N);
%nodes in the grid
for i = 2:N-1
    %down
    if mod(i-1,L) ~= 0
        num = -200;
        while sigma*num < -1
            num = randn;
        end
        Ic_zero(i,i+1) = 1 + sigma*num;
        Ic_zero(i+1,i) = Ic_zero(i,i+1);
        
        Aij(i,i+1) = 1;
        Aij(i+1,i) = 1;
    end
    %right
    if i < N - L
        num = -200;
        while sigma*num < -1
            num = randn;
        end
        Ic_zero(i,i+L) = 1 + sigma*num;
        Ic_zero(i+L,i) = Ic_zero(i,i+L);
        
        Aij(i,i+L) = 1;
        Aij(i+L,i) = 1;
    end
end

% supernodes
for i = 1:L
    num = -200;
    while sigma*num < -1
        num = randn;
    end
    Ic_zero(1,i+1) = 1 + sigma*num;
    Ic_zero(i+1,1) = Ic_zero(1,i+1);
    
    Aij(1,i+1) = 1;
    Aij(i+1,1) = 1;

    num = -200;
    while sigma*num < -1
        num = randn;
    end

    Ic_zero(N,N-i) = 1 + sigma*num;
    Ic_zero(N-i,N) = Ic_zero(N,N-i);
    
    Aij(N,N-i) = 1;
    Aij(N-i,N) = 1;
end

% for i = 2:N-1
%     %down
%     if mod(i-1,L) ~= 0
%         num = randn;
%         Ic_zero(i,i+1) = abs(1 + num);
%         Ic_zero(i+1,i) = Ic_zero(i,i+1);
%     end
%     %right
%     if i < N - L
%         num = randn;
%         Ic_zero(i,i+L) = abs(1 + num);
%         Ic_zero(i+L,i) = Ic_zero(i,i+L);
%     end
% end
% 
% % supernodes
% for i = 1:L
%     num = randn;
%     Ic_zero(1,i+1) = abs(1 + num);
%     Ic_zero(i+1,1) = Ic_zero(1,i+1);
%     
%     num = randn;
%     Ic_zero(N,N-i) = abs(1 + num);
%     Ic_zero(N-i,N) = Ic_zero(N,N-i);
% end