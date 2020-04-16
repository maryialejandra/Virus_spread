import numpy as np
import matplotlib.pyplot as mplt

N=100
L=20
ratio=0.5
iterations=200
class virus_evol:
    #ratio=
    #beta=
    #L=
    #N_part=
    def __init__(self,r0,v0,infected,time):
        self.r = r0
        self.v = v0
        self.recovery=[]
        self.infected=infected
        self.time=time
    
    def distance(self,i,j,rt):
        d=np.sqrt(sum((self.r[i]-self.r[j])**2))
        #print(d)
        return d <= rt*2
    def dislim(self,i,rt,L):
        #l1=np.absolute(r1-rt<=rt*2w_i
        l1 = (self.r[i]-rt)<=0
        #l2=np.absolute(r1-L)<=rt*2
        l2 = (self.r[i]+rt)>=L
        if l1.all():
            self.r[i] = np.ones(2)*rt
            self.v[i] = -self.v[i]
        elif l1.any():
            if l1[0]:
                self.r[i][0] = rt
                self.v[i][0] = -self.v[i][0]
                #    v1=-v1
            elif l1[1]:
                self.r[i][1] = rt
                self.v[i][1] = -self.v[i][1]
        if l2.all():
            self.r[i] = np.ones(2)*(L-rt)
            self.v[i] = -self.v[i]
        elif l2.any():
            if l2[0]:
                self.r[i][0] = L-rt
                self.v[i][0] = -self.v[i][0]
            elif l2[1]:
                self.r[i][1] = L-rt
                self.v[i][1] = -self.v[i][1]
        return self.v[i]
    
    def vel_evol(self,i,j):
        self.v[j] = -self.v[i]
        self.v[i] = -self.v[j]
        return self.v[i],self.v[j]
    def x_evol(self):
        self.r = self.v+self.r
        return self.r
    def is_infected(self,i,j):
        if (i == self.infected).any() and (j == self.infected).any():
            pass
        elif (i == self.infected).any() or (j == self.infected).any():
            if (i == self.infected).any():
                if (np.array(self.recovery) == j).any():
                    self.infected = self.infected
                else:
                    self.infected = np.append(self.infected,j)
            else:
                if (np.array(self.recovery) == i).any():
                    self.infected = self.infected          
                else:
                    self.infected = np.append(self.infected,i)
        return self.infected    
    def infected_time(self,d_infected):
        #recovery=[]
        for i in self.infected:
            #print(i)
            self.time[i] += 1
            if self.time[i] == d_infected:
                #print('se recuperó')
                self.infected = self.infected[self.infected != i]
                self.recovery.append(i)
                self.time[i] = 0
        #print('nuevos infectados',new_infected)
        return self.infected,self.recovery,self.time


count_infected=[]
count_recovery=[]
parti=[]
parti2=[]

##Funciones para inicialización
def initial_pos(N_parti,L):
    return np.random.random((N_parti,2))*L
def V0(N_p):
    return np.random.uniform(low=-0.5,high=0.5,size=(N_p,2))
def V0_0(N_p,N_v):
    vel= np.random.uniform(low=-0.5,high=0.5,size=(N_v,2))
    vel0=np.zeros((N_p-N_v,2))
    return np.concatenate((vel,vel0),axis=0)
def i_infected(N_p,N_i):
    time = np.zeros(N_p)
    infected = np.arange(0,N_i)#np.random.random_integers(low = 0,high = N_p-1,size = N_i)
    for j in infected:
        time[j] = 1
    return time, infected

r0=initial_pos(N,L)
#v0=V0_0(N,10)
#print(v0)
v0=V0(N)
time,infected=i_infected(N,1)
particulas=virus_evol(r0,v0,infected,time)
beta=0.3
d_inf=30
#print(m.r)

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
                particulas.v[i],particulas.v[j]=particulas.vel_evol(i,j)
            else:
                pass
    particulas.infected,particulas.recovery,particulas.time=particulas.infected_time(d_inf)
    if len(infected)==0:
        print('there is not spread')
        break
    count_infected.append(len(particulas.infected))
    count_recovery.append(len(np.array(particulas.recovery)))
    particulas.r=particulas.x_evol()
    #parti.append(r[0])
    #parti2.append(r[1])     



for i in range(10):
    mplt.plot(np.arange(len(count_infected)),count_infected,label='Infectados')
    mplt.plot(np.arange(len(count_recovery)),np.array(count_recovery),label='Recuperados',color='g')
    mplt.grid()
    mplt.legend()
    mplt.xlabel('Días')
    mplt.ylabel('Contagiados Activos')
    mplt.savefig('hola'+str(i)+'.png')
#mplt.show()
