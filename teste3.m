val = csvread('valores.txt');
t = csvread('ts.txt');
t = t';
m = size(val,2);
lamb1 = (2*( 1-cos(pi/(m+1))));
lamb2 = (2*( 1-cos((m*pi/((m+1))))));
for i = 1:size(val)
    for j = 1:m
        yj = j/(m+1);
        f(i,j) = exp(-lamb1*t(i))*sin(pi*(yj))+exp(-lamb2*t(i))*sin(m*pi*(yj));
    end
end
figure(1)
stem(t,val(:,1)); xlabel('variavel t'); ylabel('x(t)'); title('Elemento 1 do vetor X (Teste 3)')
hold on
plot(t,f(:,1)); 
hold off
figure(2)
stem(t,val(:,2)); xlabel('variavel t'); ylabel('x(t)'); title('Elemento 2 do vetor X (Teste 3)')
hold on
plot(t,f(:,2)); 
hold off
figure(3)
stem(t,val(:,3)); xlabel('variavel t'); ylabel('x(t)'); title('Elemento 3 do vetor X (Teste 3)')
hold on
plot(t,f(:,3)); 
hold off
figure(4)
stem(t,val(:,4)); xlabel('variavel t'); ylabel('x(t)'); title('Elemento 4 do vetor X (Teste 3)')
hold on
plot(t,f(:,4)); 
hold off
figure(5)
stem(t,val(:,5)); xlabel('variavel t'); ylabel('x(t)'); title('Elemento 5 do vetor X (Teste 3)')
hold on
plot(t,f(:,5)); 
hold off
figure(6)
stem(t,val(:,6)); xlabel('variavel t'); ylabel('x(t)'); title('Elemento 6 do vetor X (Teste 3)')
hold on
plot(t,f(:,6)); 
hold off
figure(7)
stem(t,val(:,7)); xlabel('variavel t'); ylabel('x(t)'); title('Elemento 7 do vetor X (Teste 3)')
hold on
plot(t,f(:,7)); 
hold off