function [frame] = set_frame_for_movie_two_layers(L,N,Ti,Tc,R0, state1,  state2, T_heating_vec, R_heating_vec1, R_heating_vec2, T_cooling_vec, R_cooling_vec1, R_cooling_vec2)
    figure; hold on;
    tiledlayout(1,3);
    nexttile;
    hold on;
    plot(T_heating_vec, R_heating_vec1/R0, '-ro');
    plot(T_heating_vec, R_heating_vec2/R0, '-bo');

    if ~isempty(R_cooling_vec1)
        plot(T_cooling_vec, R_cooling_vec1/R0, '-rs');
        plot(T_cooling_vec, R_cooling_vec2/R0, '-bs');
        legend('heating Top','heating Bottom','cooling Top','cooling Bottom', 'Location', 'northwest');
    else
        legend('heating Top','heating Bottom', 'Location', 'northwest');
    end
    title('blue - superconductor, green - intermediate, red - normal')
    xlim([0 Tc]);
    ylim([0,1.2]);
    ylabel('R_N(\Omega)');
    xlabel('T(K)');
    
    
    nexttile; hold on;
    for i = 2:N-1

        x = floor((i-2)/L) + 1;
        y = L - round(mod((i-2),L));
        % down
        if mod(i-1,L) ~= 0
            if state1(i,i+1) == 1
                plot([x,x],[y,y-1],'r-','LineWidth',3);
            elseif state1(i,i+1) == 2
                plot([x,x],[y,y-1],'g-','LineWidth',3);
            else
                plot([x,x],[y,y-1],'b-','LineWidth',3);
            end
        end

        %right

        if i < N - L ~= 0
            if state1(i,i+L) == 1
                plot([x,x+1],[y,y],'r-','LineWidth',3);
            elseif state1(i,i+L) == 2
                plot([x,x+1],[y,y],'g-','LineWidth',3);
            else
                plot([x,x+1],[y,y],'b-','LineWidth',3);
            end
        end
    end
    xlim([0 L+1]);
    ylim([0 L+1]);
    %title('Layer1: blue - superconductor, green - intermediate, red - normal');
    title('Top');

    
    nexttile; hold on;
    for i = 2:N-1

        x = floor((i-2)/L) + 1;
        y = L - round(mod((i-2),L));
        % down
        if mod(i-1,L) ~= 0
            if state2(i,i+1) == 1
                plot([x,x],[y,y-1],'r-','LineWidth',3);
            elseif state2(i,i+1) == 2
                plot([x,x],[y,y-1],'g-','LineWidth',3);
            else
                plot([x,x],[y,y-1],'b-','LineWidth',3);
            end
        end

        %right

        if i < N - L ~= 0
            if state2(i,i+L) == 1
                plot([x,x+1],[y,y],'r-','LineWidth',3);
            elseif state2(i,i+L) == 2
                plot([x,x+1],[y,y],'g-','LineWidth',3);
            else
                plot([x,x+1],[y,y],'b-','LineWidth',3);
            end
        end
    end
    xlim([0 L+1]);
    ylim([0 L+1]);
    %title('Layer2: blue - superconductor, green - intermediate, red - normal');
    title('Bottom');
    
    
    
    set(gcf, 'units', 'normalized');
    set(gcf, 'Position', [0, 0.1, 1, 0.5]); 
    frame = getframe(gcf);