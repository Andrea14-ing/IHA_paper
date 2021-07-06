function [o, va, vb, ma, mb]= ellipse_fit(p,beta)
% p:	[nFx2] traiettoria di cui calcolare l’ellisse di confidenza
% beta:	[scalare] livello di confidenza. es 0.95
% o:	[2x1] centro
% va, vb:	[2x1] versori di a e b
% ma, mb:	[2x1] lunghezze di a e b


% Andrea Ancillao


nF=size(p,1);          %numero di frame della traiettoria
o=mean(p,1);           %calcola il baricentro della traiettoria
p=p-repmat(o,[nF,1]);  %ridefinisce p come traiettoria rispetto al baricentro

C=(p'*p)/(nF-1);       %calcola la matrice di covarianza [2x2]
[U S V]=svd(C);        % decomposizione SVD
va=U(:,1);             % estrae l’autovettore maggiore da U
vb=U(:,2);             % estrae l’autovettore minore da U
sigma_a=sqrt(S(1,1));  % calcola sigma_a
sigma_b=sqrt(S(2,2));  % calcola sigma_b

z_beta=chi2inv(beta,2);   % distribuzione chi^2 al beta% di confidenza (2 gradi di libertà!!)
ma=sqrt(z_beta)*sigma_a; % lunghezza del semiasse maggiore
mb=sqrt(z_beta)*sigma_b; % lunghezza del semiasse minore

end

