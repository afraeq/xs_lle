function FH_Flash_pso (chi_input)

format long;

%% pra plotar maneiro

lw = 2.2;       % LineWidth
alw = 1.2;      % AxesLineWidth
msz = 8;        % MarkerSize
fsz = 13;       % Fontsize

set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultAxesLineWidth',alw);  % set the default axis line width to alw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultAxesFontSize',fsz);   % set the default font size to fsz

%% condicoes

% numero de fases
nf = 2;

% tamanhos de cadeia

n_sol = 1;
n_pol = [10 100 200 300 400 500 600 700 800 900 1000 1250 1500 1750 2000 ...
         2500 3000 4000 5000 6000 7000 8000 9000 10000];
     
n = [n_sol n_pol];

% numero de componentes
nc = length(n)

% parametro shulz-flory
p = 0.999;

% volumes molares
v_bar = n;

% volumes globais
vi_sol = 0.98;
vi_pol = (1-vi_sol)*(shulz_flory(n_pol,p)/sum(shulz_flory(n_pol,p)));
vi = [vi_sol vi_pol];

% plotando a alimentacao
plot(n_pol,vi_pol/(1-vi_sol),'--k*');

% parametros de interacao
chi = zeros(nc,nc);

for m=2:nc
    chi(1,m) = chi_input;
    chi(m,1) = chi_input;
end

%% funcao objetivo para calculo de equilibrio: energia livre de gibbs

    function deltaG = deltaG(vij_vetor)
        
        deltaG = 0;
        
        vij = zeros(nc,nf-1);
        phi_ij = vij;
        mi = vij;
        v_i_nf = zeros(1,nc);
        mi_nf = v_i_nf;
        
        for i=1:nc
            for j=1:nf-1
                vij(i,j) = vij_vetor(j+(i-1)*(nf-1));
            end
        end
        
        for j = 1:nf-1
            for i = 1:nc
                phi_ij(i,j) = vij(i,j)/sum(vij(:,j));
            end
            for i=1:nc
                sum1 = 0; sum2 = 0; sum3 = 0;
                for k=1:nc
                    if (k==i)
                        continue
                    end
                    sum1 = sum1 + (1-n(i)/n(k))*phi_ij(k,j);
                    sum2 = sum2 + phi_ij(k,j)*chi(i,k);
                    for l=1:(k-1)
                        if (l==i)
                            continue
                        end
                        sum3 = sum3 + phi_ij(k,j)*phi_ij(l,j)*chi(l,k);
                    end
                end
                mi(i,j) = log(phi_ij(i,j)) + sum1 + n(i)*((1-phi_ij(i,j))*sum2 - sum3);
                deltaG = deltaG + (vij(i,j)/v_bar(i))*mi(i,j);
            end
        end
        
        for i=1:nc
            v_i_nf(i) = vi(i) - sum(vij(i,:));
        end
        
        phi_i_nf = v_i_nf/sum(v_i_nf);
        
        for i=1:nc
            sum1 = 0; sum2 = 0; sum3 = 0;
            for k=1:nc
                if (k==i)
                    continue
                end
                sum1 = sum1 + (1-n(i)/n(k))*phi_i_nf(k);
                sum2 = sum2 + phi_i_nf(k)*chi(i,k);
                for l=1:(k-1)
                    if (l==i)
                        continue
                    end
                    sum3 = sum3 + phi_i_nf(k)*phi_i_nf(l)*chi(l,k);
                end
            end
            mi_nf(i) = log(phi_i_nf(i)) + sum1 + n(i)*(((1-phi_i_nf(i))*sum2 - sum3));
            deltaG = deltaG + (v_i_nf(i)/v_bar(i))*mi_nf(i);
        end
        
    end

%%    encontrando o minimo global usando o metodo PSO

options = optimoptions(@particleswarm,'Display','Iter','TolFun',1e-8);

bounds1 = zeros(nc,1);
bounds2 = bounds1;

for m=1:nc
    bounds1(m) = 1e-50;
    bounds2(m) = vi(m) - 1e-10;
end

[vi_fase1, deltaG_min] = particleswarm (@deltaG,(nc*(nf-1)),bounds1,bounds2,options);

%% arrumando resultados

vi
vi_fase1
length(vi)
length(vi_fase1)

vi_fase2 = (vi-vi_fase1);

phi_fase1 = vi_fase1/sum(vi_fase1)
phi_fase2 = vi_fase2/sum(vi_fase2)

%% calculando o xs e plotando

hold on

if phi_fase1(1) > phi_fase2(1)
    xs = sum(vi_fase1(2:end))/(1-vi_sol);
    plot(n_pol,phi_fase1(2:end)/sum(phi_fase1(2:end)),'--*r')
    plot(n_pol,phi_fase2(2:end)/sum(phi_fase2(2:end)),'--*b')
else
    xs = sum(vi_fase2(2:end))/(1-vi_sol);
    plot(n_pol,phi_fase2(2:end)/sum(phi_fase2(2:end)),'--*r')
    plot(n_pol,phi_fase1(2:end)/sum(phi_fase1(2:end)),'--*b')
end

legend('Alimentacao','Fase solvente','Fase polimero')
xlabel('Tamanho de cadeia')
ylabel('Fracao volumetrica')

%chi_input_str = num2str(chi_input*100,2);
%save(strcat('FH_minGibbs_pso_',chi_input_str,'.mat'),'xs','deltaG_min', 'vi_fase1', 'vi_fase2')

%%
function psi = shulz_flory (r,p)
    psi = r.*((1-p).^2).*p.^(r-1);
end

end

