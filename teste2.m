val = csvread('valores.txt');
t = csvread('ts.txt');
t = t';

for i = 1:size(val)
    f(i) = exp(-t(i))*sin(t(i))+exp(-3*t(i))*cos(3*t(i));
end
for i = 1:size(val)
    f2(i) = exp(-t(i))*cos(t(i))+exp(-3*t(i))*sin(3*t(i));
end
for i = 1:size(val)
    f3(i) = -exp(-t(i))*sin(t(i))+exp(-3*t(i))*cos(3*t(i));
end
for i = 1:size(val)
    f4(i) = -exp(-t(i))*cos(t(i))+exp(-3*t(i))*sin(3*t(i));
end
figure(1)
stem(t,val(:,1)); xlabel('variavel t'); ylabel('x(t)'); title('Elemento 1 do vetor X (Teste 2)')
hold on
plot(t,f); 
hold off

figure(2)
stem(t,val(:,2)); xlabel('variavel t'); ylabel('x(t)'); title('Elemento 2 do vetor X (Teste 2)')
hold on
plot(t,f2) 
hold off

figure(3)
stem(t,val(:,3)); xlabel('variavel t'); ylabel('x(t)'); title('Elemento 3 do vetor X (Teste 2)');
hold on
plot(t,f3)
hold off

figure(4)
stem(t,val(:,4)); xlabel('t'); ylabel('x(t)'); title('Elemento 4 do vetor X (Teste 2)');
hold on
plot(t,f4);
hold off

