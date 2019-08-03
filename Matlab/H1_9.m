function H1_9

%% dados experimentais

% massa molar media
MW = 582261;

% xs, a fracao do polimero extraida pelo xileno

xs_fraction = 3.98e-2;

% parametros das distribuicoes dos produtos finais (alimentacao do flash)

teta1_feed = 1.29e-4;
teta2_feed = 6.08e-4;
alpha_feed = .868;
q1_feed = exp(-teta1_feed);
q2_feed = exp(-teta2_feed);

% parametros das distribuicoes dos soluveis

teta1_xs = 0.000000;
teta2_xs = 6.348e-4;
alpha_xs = 0.000;
q1_xs = exp(-teta1_xs);
q2_xs = exp(-teta2_xs);

% tamanhos de cadeia

r_min = 10;
r_max = 50000;
r = 1:r_max;

% distribuicoes experimentais

feed_exp = shulz_flory2 (r, q1_feed, q2_feed, alpha_feed);
xs_exp = shulz_flory2 (r, q1_xs, q2_xs, alpha_xs);
ins_exp = (feed_exp-xs_fraction*xs_exp)/(1-xs_fraction);
ins_exp((ins_exp<0)) = 0;

% escolhendo os pseudocomponentes (classes)

n_class = 24;
delta_log_r = (log(r_max)- log(r_min))/(n_class-1);
r_class = zeros (1,n_class);
for f=1:n_class
    r_class(f) = floor(exp((f-1)*delta_log_r+log(r_min)));
end

% distribuicoes com 'n_class' componentes

feed_exp_class = zeros (1,n_class);
xs_exp_class = feed_exp_class;
ins_exp_class = feed_exp_class;

feed_exp_class(1) = sum (feed_exp(1:r_class(1)));
xs_exp_class(1) = sum (xs_exp(1:r_class(1)));
ins_exp_class(1) = sum (ins_exp(1:r_class(1)));

for f=2:n_class
    feed_exp_class(f) = sum(feed_exp(r_class(f-1):r_class(f)));
    xs_exp_class(f) = sum(xs_exp(r_class(f-1):r_class(f)));
    ins_exp_class(f) = sum(ins_exp(r_class(f-1):r_class(f)));
end

feed_exp_class = feed_exp_class/sum(feed_exp_class);
xs_exp_class = xs_exp_class/sum(xs_exp_class);
ins_exp_class = ins_exp_class/sum(ins_exp_class);

% transladando as distribuicoes

aux = r_class;

aux(1) = floor(exp(((log(r_class(1))+log(1))/2)));

for f=2:n_class
    aux(f)=floor(exp(((log(r_class(f))+log(r_class(f-1)))/2)));
end

r_class = aux;

%% condicoes

% numero de fases

nf = 2;

% tamanhos de cadeia

n_sol = 1;
n_pol = r_class;
n = [n_sol n_pol];

% numero de componentes

nc = length(n);

% volumes molares

v_bar = n;

% volumes globais

vi_sol = 0.97738;
vi_pol = (1-vi_sol)*feed_exp_class;
vi = [vi_sol vi_pol];

%% rodando a estimacao de parametros

% inicializando o chi

chi = zeros(nc,nc);

% faixa das variaveis de otimizacao

bounds1_eq = zeros(nc,1);
bounds2_eq = bounds1_eq;

for m=1:nc
    bounds1_eq(m) = 1e-50;
    bounds2_eq(m) = vi(m) - 1e-10;
end

% estima!

options_eq = optimoptions(@particleswarm,'Display','Off','TolFun',1e-8,'UseParallel',true);
options_estima = optimoptions(@particleswarm,'Display','Off','TolFun',1e-4);
chi_otimo = particleswarm (@objF,1,0.50,0.55,options_estima);

% calculando o equilibrio para o chi otimo encontrado

for m=2:nc
    chi(1,m) = chi_otimo;
    chi(m,1) = chi_otimo;
end

% composicoes otimas

vi_fase1_otimo = particleswarm (@deltaG,(nc*(nf-1)),bounds1_eq,bounds2_eq,options_eq);
vi_fase2_otimo = (vi-vi_fase1_otimo);

phi_fase1_otimo = vi_fase1_otimo/sum(vi_fase1_otimo);
phi_fase2_otimo = vi_fase2_otimo/sum(vi_fase2_otimo);

% salvando arquivo com os dados e resultados a serem processados posteriormente

save(file_name, 'xs_fraction', 'MW', 'chi_otimo', 'phi_fase1_otimo', 'phi_fase2_otimo','n_pol','feed_exp_class','xs_exp_class','ins_exp_class');

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

%% funcao objetivo para a estimacao de parametros

    function objF = objF(input)
        
        for o=2:nc
            chi(1,o) = input;
            chi(o,1) = input;
        end
        
        % encontrando o minimo global usando o metodo PSO
        
        vi_fase1 = particleswarm (@deltaG,(nc*(nf-1)),bounds1_eq,bounds2_eq,options_eq);
        
        % arrumando resultados
        
        vi_fase2 = (vi-vi_fase1);
        
        phi_fase1 = vi_fase1/sum(vi_fase1);
        phi_fase2 = vi_fase2/sum(vi_fase2);

	% calculando a funcao objetivo
        
        if (phi_fase1(1)>phi_fase2(1))
            xs_modelo = vi_fase1(2:end)/sum(vi_fase1(2:end));
            ins_modelo = vi_fase2(2:end)/sum(vi_fase2(2:end));
        else
            ins_modelo = vi_fase1(2:end)/sum(vi_fase1(2:end));
            xs_modelo = vi_fase2(2:end)/sum(vi_fase2(2:end));
        end
        
        objF = sum((xs_modelo-xs_exp_class).^2)+sum((ins_modelo-ins_exp_class).^2);
        
    end

%% distribuicao de shulz-flory

    function psi = shulz_flory2 (r,q1,q2,alpha)
        psi = (alpha*(1-q1)*(1-q1)*(q1.^(r-1))).*r + ((1-alpha)*(1-q2)*(1-q2)*(q2.^(r-1))).*r; 
    end

end
