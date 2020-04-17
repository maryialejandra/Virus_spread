import numpy as np
import matplotlib.pyplot as mplt

N=200 # Numero de personas
L=20  #Logitud de la caja
ratio=0.5 #Radio de las particulas
iterations=200 #Iteraciones del código 
class virus_evol: 

    def __init__(self,r0,v0,infected,time): #Inicialización de las variables
        self.r = r0
        self.v = v0
        self.recovery=[]
        self.infected=infected
        self.time=time
    
    def distance(self,i,j,rt): # Mide la distancia entre particulas
        d=np.sqrt(sum((self.r[i]-self.r[j])**2))
        return d <= rt*2
    def dislim(self,i,rt,L):  #Modificación de velocidad cuando la particula choca contra las paredes
        l1 = (self.r[i]-rt)<=0
        l2 = (self.r[i]+rt)>=L
        if l1.all(): #Choque en las esquinas
            self.r[i] = np.ones(2)*rt
            self.v[i] = -self.v[i]
        elif l1.any():
            if l1[0]: #Choque en el eje x
                self.r[i][0] = rt
                self.v[i][0] = -self.v[i][0]
                #    v1=-v1
            elif l1[1]: #Choque en el eje y
                self.r[i][1] = rt
                self.v[i][1] = -self.v[i][1]
        if l2.all():  #Choque en las esquinas
            self.r[i] = np.ones(2)*(L-rt)
            self.v[i] = -self.v[i]
        elif l2.any():
            if l2[0]:  #Choque en el eje x
                self.r[i][0] = L-rt
                self.v[i][0] = -self.v[i][0]
            elif l2[1]:  #Choque en el eje y
                self.r[i][1] = L-rt
                self.v[i][1] = -self.v[i][1]
        return self.v[i]

    def vel_evol(self,i,j): #Velocidad después de la colisión elástica
        self.v[j] = -self.v[i]
        self.v[i] = -self.v[j]
        return self.v[i],self.v[j]
    def x_evol(self): #Evolución de las posiciones
        self.r = self.v+self.r
        return self.r
    def is_infected(self,i,j): #Define si las particulas se infectan y cuales
        if (i == self.infected).any() and (j == self.infected).any(): #Significa que ambas particulas estan infectadas
            pass
        elif (i == self.infected).any() or (j == self.infected).any(): 
            if (i == self.infected).any(): #La particula i está infectada entonces infectaría a la particula j
                if (np.array(self.recovery) == j).any(): #Si j está recuperada no la infecta
                    self.infected = self.infected
                else:
                    self.infected = np.append(self.infected,j)
            else: #La particula j está infectada entonces infectaría a la particula i
                if (np.array(self.recovery) == i).any(): 
                    self.infected = self.infected   #Si i está recuperada no la infecta       
                else:
                    self.infected = np.append(self.infected,i)
        return self.infected    
    def infected_time(self,d_infected): #Define cuanto tiempo lleva alguien enfermo
        for i in self.infected:
            self.time[i] += 1
            if self.time[i] == d_infected: #Si llega al tiempo de recuperación, cambia su estado a recuperada
                self.infected = self.infected[self.infected != i]
                self.recovery.append(i)
                self.time[i] = 0
        return self.infected,self.recovery,self.time



np.random.seed(6) #Para realizar nuestras comparaciones con las mismas condiciones iniciales
##Funciones para inicialización de condiciones
def initial_pos(N_parti,L):
    return np.random.random((N_parti,2))*L
def V0(N_p): # Velocidad sin medidas de aislamiento 
    return np.random.uniform(low=-0.5,high=0.5,size=(N_p,2))
def V0_0(N_p,N_v): # Velocidad con medidas de aislamiento 
    vel= np.random.uniform(low=-0.5,high=0.5,size=(N_v,2))
    vel0=np.zeros((N_p-N_v,2))
    return np.concatenate((vel,vel0),axis=0)
def i_infected(N_p,N_i): #Inicialización de infectados
    time = np.zeros(N_p)
    infected = np.arange(0,N_i)
    for j in infected:
        time[j] = 1
    return time, infected

#r0=initial_pos(N,L)
#v0=V0_0(N,50)
#v0=V0(N)
#time,infected=i_infected(N,1)
#particulas=virus_evol(r0,v0,infected,time)
beta=0.7 # Probabilidad de contagio 
d_inf=21
betas=np.array([0.3,0.4,0.5,0.7])

def evolution(particulas,beta):
    count_infected=[]
    count_recovery=[]
    for m in range(iterations):
        for i in range(N):
            particulas.v[i]=particulas.dislim(i,ratio,L)
            for j in range(N):
                if i <= j:  
                    break
                l=particulas.distance(i,j,ratio)
                if l:
                    if np.random.random(1)<=beta:
                        particulas.infected=particulas.is_infected(i,j)
                    #particulas.v[i],particulas.v[j]=particulas.vel_evol(i,j) # Esta se comenta  para no tener cambio de velocidades en un choque
                else:
                    pass
        particulas.infected,particulas.recovery,particulas.time=particulas.infected_time(d_inf)
        if len(infected)==0:
            print('there is not spread')
            break
        count_infected.append(len(particulas.infected))
        count_recovery.append(len(np.array(particulas.recovery)))
        particulas.r=particulas.x_evol()
    return particulas.r,count_infected,count_recovery

for N2 in range(25,150,25): # Se realiza la gráfica para observar el comportamiento de la medida de aislamiento.
    r0=initial_pos(N,L)
    v0=V0_0(N,N2)

    time,infected=i_infected(N,1)

    particulas=virus_evol(r0,v0,infected,time)
    _,count_infected,count_recovery=evolution(particulas,beta)
    mplt.plot(np.arange(len(count_infected)),count_infected,label='Personas en movimiento='+str(N2))
#mplt.plot(np.arange(len(count_recovery)),np.array(count_recovery),label='Recuperados',color='g')
mplt.grid()
mplt.legend()
mplt.xlabel('Días')
mplt.ylabel('Contagiados Activos')
mplt.xlim(0,200)
mplt.ylim(0,200)
mplt.savefig('fig8.png')
mplt.show()


dias = count_infected.index(max(count_infected))
maximo = max(count_infected)
print(dias,maximo)


minimo = min(count_infected)
diasmin = count_infected.index(min(count_infected))

print(diasmin,minimo)                                                                                                                                                                                                                                                                                 
#mplt.savefig('hola'+str(i)+'.png')
#mplt.show()
