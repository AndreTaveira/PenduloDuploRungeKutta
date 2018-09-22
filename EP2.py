import math

global var_rebelde
var_rebelde = 7
global teste1
teste1 = False
global teste2
teste2 = False
global teste3
teste3 = False

def main():

    global teste1
    teste1 = False
    global teste2
    teste2 = False
    global teste3
    teste3 = False
    
    ''' decide se cria arquivo txt para teste em matlab '''
    imprime = True
    
    ''' seleciona se deseja realizar os testes propostos ou o teste do pendulo'''
    selecao = int(input("digite: \n * '1' para realizar um dos testes sugeridos no enunciado: \n * '2' para realizar teste do pêndulo duplo: \n   "))
    print("")
    # observe que as variaveis xi ti h de retorno da função runge kutta não são usadas aqui.
    # Foram criadas apenas para testes durante a realização deste EP
    if selecao == 1:
        teste = int(input("digite o teste que deseja realizar ( '1' , '2' ou '3' ) : "))
        print("")
        if teste == 1:
            teste1 = True
            xi,ti,h,val,ts = RKF45([-18.95],1.05,0.1,1e-5,f_teste1,3)

        elif teste == 2:
            teste2 = True
            xi,ti,h,val,ts = RKF45([1,1,1,-1],0,0.1,1e-5,f_teste2,2)

        elif teste == 3:
            teste3 = True
            var_rebelde = int(input("Digite o tamanho m da matriz de teste (tamanho 7 sugerido) : "))
            print("")
            x0 = []
            ''' cria x0 inicial para este teste'''
            for z in range(1,var_rebelde+1,1):
                yz = z/(var_rebelde+1)
                x0.append( math.sin(math.pi*yz) + math.sin(var_rebelde*math.pi*yz) )
            xi,ti,h,val,ts = RKF45(x0,0,0.1,1e-5,f_teste3,2)
    # teste do pendulo 
    if selecao == 2:
        th1 = float(input("Digite Theta 1 : "))
        th2 = float(input("Digite Theta 2 : "))
        w1 = 0
        w2 = 0
        conv = math.pi/180
        x0 = [th1*conv,th2*conv,w1*conv,w2*conv]
        xi,ti,h,val,ts = RKF45(x0,0,0.1,1e-8,f_tarefa,120)
        # calcula agora a energia do sistema pra cada elemento dos valores calculados
        T = []
        V = []
        E = []
        for i in range(0,len(val),1):
            T.append((val[i][2]**2) + ((math.cos(val[i][0]-val[i][1]))*(val[i][2])*(val[i][3]))+((val[i][3])**2)/2)
            V.append(2*(1-math.cos(val[i][0])) + (1-math.cos(val[i][1])) )
            E.append(T[i]+V[i])

        arqE = open("energia.txt",'w')
        for i in range(len(E)):
            arqE.write( str ( E[i] ) )
            arqE.write('\n')
        arqE.close()
        L = 1
        plot_xy = []
        # cria arquivo que armazena as posições (x1,y1) e (x2,y2) 
        for i in range (0,len(val)):
            pto = []
            pto.append(L*math.sin(val[i][0]))
            pto.append(-L*math.cos(val[i][0]))
            pto.append(L*math.sin(val[i][0])+L*math.sin(val[i][1]))
            pto.append(-L*math.cos(val[i][0])-L*math.cos(val[i][1]))
            plot_xy.append(pto)
        escreve_pt(plot_xy)
        
    print ("Valor de x(%.2f) calculado: " %(ti), xi , end="\n")

    if imprime == True:
        escreve(val,ts)
    
def f_teste1(t,x):
    
    return [( 1 + (x[0]-t)**2 )]

def f_teste2 (t,x):

    ''' estamos assumindo que x é um vetor de tamanho 4x1 '''
    A = [[-2,-1,-1,-2],[1,-2,2,-1],[-1,-2,-2,-1],[2,-1,1,-2]]
    f_ret = [0,0,0,0]
    
    for i in range (len(x)):
        for j in range(len(x)):
            f_ret[i] += A[i][j]*x[j]     

    return f_ret

def f_teste3(t,x):

    ''' Note que estamos realizando a multiplicação de uma matriz
    por um vetor tridiagonal de for eficiente, ignorando onde há zeros '''
    vec_ret = []
    for i in range (var_rebelde):
        vec_ret.append(0)
        for j in range (max(0,i-1),min(var_rebelde-1,i+1)+1,1):
            if j == i-1:
                vec_ret[i] += x[j]
            if j == i:
                vec_ret[i] += -2*x[j]
            if j == i+1:
                vec_ret[i] += x[j]
                
    return vec_ret

def f_tarefa(t,x):

    ''' realiza a aplicação do sistema que descreve o pendulo duplo '''
    vec_ret = [0,0,0,0]
    vec_ret[0] = x[2]
    vec_ret[1] = x[3]
    vec_ret[2] = ( math.cos(x[0]-x[1])*math.sin(x[1]) - 2*math.sin(x[0]) - math.sin(x[0]-x[1])*(math.cos(x[0]-x[1])*((x[2])**2) + ((x[3])**2)))/(1+(math.sin(x[0]-x[1]))**2)
    vec_ret[3] = ( 2*(math.cos(x[0]-x[1])*math.sin(x[0])-math.sin(x[1])) + math.sin(x[0]-x[1])*(math.cos(x[0]-x[1])*((x[3])**2) + 2*((x[2])**2)))/(1+(math.sin(x[0]-x[1]))**2)
    return vec_ret


def atualiza_h (alpha,h,tf,ti):
    ''' Aqui atualizamos o valor de h para 0.001 < h < 1 e igual ao
    complemento a tf (tempo final) caso esteja muito proximo'''
    h_linha = alpha*h
    h = min(h_linha,1.0)
    h = max(h_linha,0.001)
    h = min(h,tf-ti)
    
    
    return h

def calcula_alpha(eps,h,xtil,x):
    ''' realiza o calculo de alpha. Note que foi utilizado um valor minimo no denominador para evitar que ocorra divisãopor zero '''
    den = max(( 2*max ( modulo ( vec_sum ( xtil , ppEsc(x,-1) ) ) ) ),0.0000000000000000001)
    alpha = pow ( (eps * h / den ) , 1/4 )

    return alpha
    

def kis(ti,xi,h,f):
    ''' Calculo dos k's com auxilio das funções ppEsc que faz o produto de um
    escalar por um vetor e a função vec_sum que soma dois vetores'''
    k1 = ppEsc(f( ti , xi ),h)
    k2 = ppEsc(f( ti + h/4 , vec_sum ( xi , ppEsc( k1 , 1/4 ) ) ),h)
    k3 = ppEsc(f( ti + 3*h/8 , vec_sum ( xi , vec_sum( ppEsc(k1,3/32) , ppEsc(k2,9/32) ) ) ),h)
    k4 = ppEsc(f( ti + 12*h/13 , vec_sum ( xi , vec_sum( ppEsc( k1 , 1932/2197 ) , vec_sum( ppEsc( k2 , -7200/2197 ) , ppEsc( k3 , 7296/2197 ) ) ) )),h)
    k5 = ppEsc(f( ti + h , vec_sum ( xi , vec_sum(ppEsc(k1,439/216) , vec_sum( ppEsc(k2,-8) , vec_sum(ppEsc(k3,3680/513) , ppEsc(k4,-845/4104)))))),h)
    k6 = ppEsc(f( ti + h/2 , vec_sum ( xi , vec_sum (ppEsc(k1,-8/27) , vec_sum( ppEsc(k2,2) , vec_sum( ppEsc(k3,-3544/2565) , vec_sum( ppEsc(k4,1859/4104) , ppEsc(k5,-11/40) ) ) ) ) )),h)

    ki = [k1,k2,k3,k4,k5,k6]

    return ki

def vec_sum(a,b):
    ''' Faz a soma de dois vetores '''
    #supoe vetores de mesmo tamanho! Usar com cuidado
    vec = []
    for i in range(len(a)):
        vec.append(a[i]+b[i])
    return vec

def ppEsc(a,k):
    ''' Faz o produto de um escalar por um vetor '''
    vec = []
    for i in range(len(a)):
        vec.append(k*a[i])
        
    return vec
    
def modulo(a):
    ''' Faz o módulo de cada elemento de um vetor '''
    vec = []
    for i in range(len(a)):
        vec.append(abs(a[i]))
    return vec

def print_matrix(matriz):
    ''' imprime um vetor '''
    print("[", end=" ")
    for n in range (0,len(matriz),1):
        print (" %f" %(matriz[n]), end="")
    print("]")
    return

def erro_t3( x , t):
    ''' função implementada para avaliar Erro dado por max ||Xi − Xi_exato|| como pedido no enunciado '''
    Xi_exat = []
    m = var_rebelde
    lamb1 = 2*(1 - math.cos(math.pi/(m + 1)))
    lamb2 = 2*(1 - math.cos(m*math.pi/(m + 1)))
    for i in range (1,m+1,1):
        yi = i/(m+1)
        Xi_exat.append(math.exp(-lamb1*t)*math.sin(math.pi*yi) + math.exp(-lamb2*t)*math.sin(m*math.pi*yi))
    erro = max ( modulo ( vec_sum ( Xi_exat, ppEsc(x,-1) ) ) )
    return erro

    
def RKF45 (x0,t0,h0,eps,f,tf):
    
    xi = x0
    ti = t0
    h = h0
    val = []
    ts = []
    hs = []
    # executa os calculos de x enquanto ti não chega ao fim 
    while( ti < tf ):
        ki = kis(ti,xi,h,f)
        # calcula x_i+1 e xtil_1+1 e o tau que avalia o erro 
        x_imais1 = vec_sum(xi , vec_sum(ppEsc(ki[0],25/216) , vec_sum(ppEsc(ki[2],1408/2565) , vec_sum(ppEsc(ki[3],2197/4104) , ppEsc(ki[4],-1/5)))))
        xtil_imais1 = vec_sum(xi , vec_sum(ppEsc(ki[0],16/135) , vec_sum(ppEsc(ki[2],6656/12825) , vec_sum(ppEsc(ki[3],28561/56430) , vec_sum(ppEsc(ki[4],-9/50) , ppEsc(ki[5],2/55))))))
        tau_imais1 = max(modulo(vec_sum(xtil_imais1 , ppEsc(x_imais1,-1))))/h
        # caso o erro esteja dentro de um valor aceitável 
        if tau_imais1 <= eps:

            # atualiza os valores para a continuação do método 
            for p in range (len(xi)):
                xi[p] = x_imais1[p]
            ti = ti + h
            alpha = calcula_alpha(eps,h,xtil_imais1,x_imais1)
            h = atualiza_h(alpha,h,tf,ti)
            val.append(x_imais1)
            ts.append(ti)
            hs.append(h)

            # as variaveis teste são globais e foram feitas para lidar com os prints exigidos nos testes do enunciado 
            if teste1 == True:
                print ("Valor de proximo h:",h, end="\n")
                print ("Valor de ti empregado:",ti, end="\n")
                print ("Valor de xi(%.3f): "%(ti),x_imais1, end="\n")
                print ("Erro entre x(t) estimado e a solução exata: %f" %( (abs(x_imais1[0] - (min(ti + h,tf) + 1/(1 - min(ti + h,tf))) )) ), end="\n\n")

            elif teste2 == True:
                print ("Valor de proximo h:",h, end="\n")
                print ("Valor de ti empregado:",ti, end="\n")
                print ("Valor de xi(%.3f): "%(ti),x_imais1, end="\n")
                print ("Erro estimado por max ||Xi − X_barra(ti)||: %f" %( max ( modulo ( vec_sum ( xtil_imais1 , ppEsc(x_imais1,-1) ) ) ) ), end="\n\n")

            elif teste3 == True:
                print ("Valor de proximo h:",h, end="\n")
                print ("Valor de ti empregado:",ti, end="\n")
                print ("Valor de xi(%.3f): "%(ti),x_imais1, end="\n")
                print ("Erro estimado por max ||Xi − X_barra(ti)||: %f" %( max ( modulo ( vec_sum ( xtil_imais1 , ppEsc(x_imais1,-1) ) ) ) ), end="\n")
                print ("Erro dado por max ||Xi − Xi_exato||: %f" %(erro_t3( x_imais1 , min(ti + h,tf))) , end="\n\n")
                
        # Caso a precisao não tenha sido atingida, atualiza o alpha e o h e refaz o teste 
        else:
            alpha = calcula_alpha(eps,h,xtil_imais1,x_imais1)
            h = atualiza_h(alpha,h,tf,ti)

    return (xi,ti,hs,val,ts)

def escreve(val,ts):
    ''' escreve em arquivos .txt os resultador obtidos com o metodo '''
    arq = open("valores.txt",'w')
    for i in range(len(val)):
        for j in range(len(val[0])):
            arq.write( str ( val[i][j] ) )
            arq.write(' ')
        arq.write('\n')
        
    arq.close()
    arq2 = open("ts.txt",'w')
    for i in range(len(ts)):
        arq2.write( str ( ts[i] ) )
        arq2.write('\n')
    arq2.close()

def escreve_pt(pontos):
    ''' função exclusiva do tarefa do pendulo para criar arquivo
    com as posições do pendulo ao longo do tempo '''
    arq = open("pontos.txt",'w')
    for i in range(len(pontos)):
        for j in range(len(pontos[0])):
            arq.write( str (pontos[i][j] ) )
            arq.write(' ')
        arq.write('\n')

main()
