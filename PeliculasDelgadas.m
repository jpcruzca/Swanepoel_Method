classdef PeliculasDelgadas
    methods
        function s = n_subst(obj,Ts)
               %{
            método que calcula el indice de refracción del sustrato desnudo a
            partir del espectro de transmitancia
               %}
            s=1./Ts+(1./Ts.^2-1).^(0.5);
        end 

        function n = n_transp(obj,s,Tm)
               %{
            método que calcula el indice de refracción de la película a
            partir del espectro de transmitancia en la región de transparencia
               %}
            M = 2*s./Tm - (s.^2+1)./2;
            n = (M+(M.^2-s.^2).^(0.5)).^(0.5);
        end 

        function d = espesor(obj,lamda,n)
               %{
        calcula el espesor del recubrimiento utilizando dos máximos consecutivos
               %}
            d = [];
            for i = 1:length(lamda)-1
                d = [d; -lamda(i+1)*lamda(i)/(2*(lamda(i)*n(i+1)-lamda(i+1)*n(i)))/2];
            end
        end

        function avg = average(obj,x)
            %promedio
            avg = sum(x)/numel(x);
        end

        function m = ordenm(obj,avgd,n,lamda)
            % calculo del orden m, ya aproximación al (semi) entero mas cercano
            m = 2*n*avgd./lamda;
            maux = round(m*2); 
            m = maux/2;
        end 

        function d_m = mespesor(obj,m,lamda,n)
            % espesor considerando el orden m
            d_m = m.*lamda./(2*n);
        end 

        function n = n_baja_media(obj,s,TM,Tm)
            %% calcula el indice de refraccion en la region baja media
            N = 2*s.*(TM-Tm)./(Tm.*TM)+(s.^2+1)/2;
            n = (N+(N.^2-s.^2).^(0.5)).^(0.5);
        end 

        function x = xTM(obj,n,s,TM)
            EM = (8.*n.^2.*s)./TM+(n.^2-1).*(n.^2-s.^2);
            x = (EM-(EM.^2-(n.^2-1).^3.*(n.^2-s.^4)).^(0.5))./((n-1).^3.*(n-s.^2));
        end 
        
        function x = xTm(obj,n,s,Tm)
            Em = (8.*n.^(2).*s)./Tm -(n.^2-1).*(n.^2-s.^2);
            x = (Em-(Em.^2-(n.^2-1).^3.*(n.^2-s.^4)).^(0.5))./((n-1).^3.*(n-s.^2));
        end 
        
        function alpha = alpha_x(obj,x,d)
            alpha = (-1/d)*log(x);
        end
        
        function x = xTi(obj,n,s,Tm,TM)
            Ti = (2*Tm.*TM)./(TM+Tm);
            F = (8*n.^2.*s)./Ti;
            x = (F-(F.^2-(n.^2-1).^3.*(n.^2-s.^4)).^(0.5))./((n-1).^3.*(n-s.^2));
        
        end
        
        function E = w_to_energy(obj,w)
            % convierte de longitud de onda nm a energía en eV
            E = 1.239841*10^3./w;
        end
        
        function tauc(obj,E,alpha)
            y = alpha.*E;
            figure
            hold on
            tipo = 2;
            y = y.^tipo*10000;
            plot(E,y,'bo');
            

            tipo = 2/3;
            y = alpha.*E;
            y = y.^tipo;
            p=polyfit(E(1:4),y(1:4),1);
            y2 = polyval(p,E(1:7));
            Eg = -p(2)/p(1)
            plot(E,y,'go')
            plot(E(1:7),y2,'r-');
            xlabel('Energía, [eV]')
            ylabel('\alpha h \nu, [eV/nm]')
            legend({'Gap directo, n=2','Gap directo, n=2/3'})
            ylim([-0.005 0.017])
            hold off
            ax = gca;
            exportgraphics(ax,'taucdir.jpg','Resolution',300)
        end 

        function tauc_ind(obj,E,alpha)
            y = alpha.*E;
            figure
            hold on
            tipo = 1/2;
            y = y.^tipo;
            plot(E,y,'bo')
            ylabel('(\alpha h \nu)^n, [eV/nm]')
            xlabel('Energía, [eV]')
            tipo = 1/3;
            y = alpha.*E;
            y = y.^tipo;
            p=polyfit(E(1:7),y(1:7),1);
            y2 = polyval(p,E(1:9));
            Eg = -p(2)/p(1)
     
            plot(E,y,'go')
            plot(E(1:9),y2,'r-');
            xlabel('Energía, [eV]')
            ylabel('(\alpha h \nu)^n, [eV/nm]')
            legend({'Gap indirecto, n=1/2','Gap Indirecto, n=1/3'})
            hold off
            ax = gca;
            exportgraphics(ax,'taucindir.jpg','Resolution',300)
        end 
        
        function k = extincion(obj,lamda,alpha)
            k = alpha.*lamda/4*pi;
        end 



        function [T w] = T_Teo(obj,lamda,n,s,k,d,x)
            w = 550:5:1100;
            pn=polyfit(1./lamda.^2,n,4);
            n = polyval(pn,1./w.^2);
            ps = polyfit(lamda,s,4);
            s = polyval(ps,w);
            pk = polyfit(1./lamda,k,6);
            k = polyval(pk,1./w);
            px = polyfit(lamda,x,6);
            x = polyval(px,w);
            lamda = w;

            phi = 4*pi*n.*d./lamda;
            A = 16*s.*(n.^2+k.^2);
            B = ((n+1).^2+k.^2).*((n+1).*(n+s.^2)+k.^2);
            C = ((n.^2-1+k.^2).*(n.^2-s.^2+k.^2)-2.*k.^2.*(s.^2+1))*2.*cos(phi)-k.*(2.*(n.^2-s.^2+k.^2)+(s.^2+1).*(n.^2-1+k.^2)).*2.*sin(phi);
            D = ((n-1).^2+k.^2).*((n-1).*(n-s.^2)+k.^2);

            T = A.*x./(B-C.*x+D.*x.^2);
        end 


    end 
end