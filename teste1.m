val = csvread('valores.txt');
t = csvread('ts.txt');
t = t';
for i = 1:size(val)
    f(i) = t(i) + 1/(1 - t(i));
end
figure(1)
stem(t,val); xlabel('variavel t'); ylabel('x(t)'); title('vetor X (Teste 1)')
hold on
plot(t,f); 
hold off
