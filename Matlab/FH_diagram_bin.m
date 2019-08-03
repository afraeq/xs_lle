function FH_diagram_bin (r1,r2)   % r1 > r2

    % ponto critico
    
    phi_C = 1/(1+(r1/r2)^0.5);
    chi_C = 0.5*(1/(r1^.5)+1/(r2^.5))^2;
    T_C = 0.5/chi_C;
    
    % vetor de composicoes da fase 1
    
    phi11_assumed = (0.01*phi_C) : (0.005*phi_C) : (phi_C-0.01*phi_C);
    
    % estimativa inicial
    
    guess = .99;
    
    % calculando a binodal
    
    for i= 1:length(phi11_assumed)
        
            phi11 = phi11_assumed(i);
            
            % encontrando as composicoes da fase 2
            
            [result,~,exitflag] = fzero(@binodal,guess,optimset('Display','Off'));

            % se o metodo numerico encontrar uma raiz complexa, mudar estimativa inicial
            
            while exitflag == -4
                guess = guess - 0.005;     
                [result,~,exitflag] = fzero(@binodal,guess);               
            end

            % composicoes da fase 2 e temperatura normalizada
            
            phi12_solved(i) = result;
            T(i) = 0.5/chi;
            
    end
    
    % calculando a espinodal
    
    x_spinodal = [phi11_assumed phi12_solved];
    T_spinodal = ((0.5./(1./(2*r1.*x_spinodal)+1./(2*r2.*(1-x_spinodal)))));
    
    % tracando a curva de deltaG
    
    subplot(2,1,1)
      
    T_tan = median(T(find(T<median(T))));   
    phi11_tan = median(phi11_assumed(find(phi11_assumed<median(phi11_assumed))));
    phi12_tan = median(phi12_solved(find(phi12_solved>median(phi12_solved))));
    
    chi_tan = 0.5/T_tan;
        
    phi_deltaG_dom = (0.01*phi_C):.001:phi12_solved(1);
        
    plot(phi_deltaG_dom, deltaG(phi_deltaG_dom),'k')
    
    YL = get(gca,'YLim');
    line([phi11_tan phi11_tan],YL,'Color','g');
    line([phi12_tan phi12_tan],YL,'Color','g');
    
    line([phi11_tan phi12_tan],[deltaG(phi11_tan) deltaG(phi12_tan)],'Color','m');
    
    xlabel('\Phi_0');
    ylabel('\Delta G');
    
    % tracando binodal e espinodal
    
    subplot(2,1,2)
    
    bin1 = plot(phi11_assumed, T,'k');
    hold on
    plot(phi12_solved, T,'k');
    crit = plot(phi_C, T_C,'*r');
    spin = plot(x_spinodal, T_spinodal,'ob','MarkerSize',1.5);
    xlabel('\Phi_0')
    ylabel('T/\theta')
    legend([bin1,spin,crit],'binodal','spinodal','critical point')
    
    XL = get(gca,'XLim');
    line(XL,[T_tan T_tan],'Color','m');
    
    YL = get(gca,'YLim');
    line([phi11_tan phi11_tan],YL,'Color','g');
    line([phi12_tan phi12_tan],YL,'Color','g');
   
    % funcao objetivo para calcular a binodal
    
    function obj = binodal (phi12)
        phi21 = 1 - phi11;
        phi22 = 1 - phi12;
        chi = (log(phi11) + (1-r1/r2)*phi21 - (1-r1/r2)*phi22 -log(phi12))/(r1*phi22.^2-r1*phi21.^2);
        obj = log(phi21) + (1-r2/r1).*(phi11-phi12) -log(phi22) +r2.*chi.*(phi11.^2-phi12.^2);
    end

    % funcao para plotar deltaG
    
    function result = deltaG (phi1)
        phi2 = 1 - phi1;
        result = ((phi1/r1).*log(phi1) + (phi2/r2).*log(phi2)) + chi_tan.*phi1.*phi2;        
    end

end
