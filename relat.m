val = csvread('valores.txt');
t = csvread('ts.txt');
E = csvread('energia.txt');
t = t';
figure(1)
plot(t,val(:,1));  xlabel('variavel t'); ylabel('x(t)'); title('vetor X1 (Teste 1)')
figure(2)
plot(t,val(:,2)); xlabel('variavel t'); ylabel('x(t)'); title('vetor X2 (Teste 1)')
figure(3)
plot(t,val(:,3)); xlabel('variavel t'); ylabel('x(t)'); title('vetor X3 (Teste 1)')
figure(4)
plot(t,val(:,4)); xlabel('variavel t'); ylabel('x(t)'); title('vetor X4 (Teste 1)')

ym = abs((E(1) - E(size(E,1))))*100
m = mean(E)
figure(5)
plot(t,E); xlabel('variavel t'); ylabel('Energia'); title('Energia com relação ao tempo')
xlim([0 120])
ylim([(m-ym) (m+ym)])

pto = csvread('pontos.txt');
figure(6)
plot(pto(:,1),pto(:,2),'r'); xlabel('x'); ylabel('y'); title(' Pontos ocupados pela posição dos pêndulos ao longo do tempo');
hold on
plot(pto(:,3),pto(:,4),'b');
hold off

figure(7)
plot3(pto(:,1),pto(:,2),t)
figure(8)
plot3(pto(:,3),pto(:,4),t)