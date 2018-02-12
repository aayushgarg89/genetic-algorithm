from random import *
from operator import *
from math import *

def individual_x1(mi,ma):
    #x1=-1.0
    x1=uniform(mi,ma)
    return x1

def individual_x2(mi,ma):
    #x2=0.003
    x2=uniform(mi,ma)
    return x2

def individual_x3(mi,ma):
    #x3=0.00006
    x3=uniform(mi,ma)
    return x3

def individual_x4(mi,ma):
    #x4=-0.00006
    x4=uniform(mi,ma)
    return x4

def individual_lambda(mi,ma):
    #lambd=17
    lambd=uniform(mi,ma)
    return lambd

def individual_rc(mi,ma):
    #rc=0.0005
    rc=uniform(mi,ma)
    return rc

def individual_b(mi,ma):
    #b=0.1
    b=uniform(mi,ma)
    return b

def fitness(x1,x2,x3,x4,lambd,rc,b):
    t,rha,rhc,pa,pc,a,ilimitden,l,ns=353.15,1,1,3,5,27,0.860,0.127,24
    sum=0
    for x in range(0,15):
        i=1.1+x
        iden=i/a
        psath2o=10**((0.0295*(t-273.15))-(0.0000919*((t-273.15)**2))+(0.000000144*((t-273.15)**3))-2.18)
        ph2=0.5*rha*psath2o*((1/((rha*psath2o/pa)*exp((1.635*iden)/t**1.334)))-1)
        #print(ph2)
        po2=rhc*psath2o*((1/((rhc*psath2o/pc)*exp((4.192*iden)/t**1.334)))-1)
        #print(po2)
        enernst=1.229-(0.85*0.001*(t-298.15))+(0.000043085*t*log1p(ph2*sqrt(po2)))
        #print(enernst)
        co2=po2/(5080000*exp(-498/t))
        nact_model=-(x1+(x2*t)+(x3*t*log1p(co2))+(x4*t*log1p(i)))
        #print(x1,(0.003*t),x3*t*log1p(co2),x4*t*log1p(i))
        nconc_model=-b*log1p(1-(iden/ilimitden))
        #print(b)
        pm_model=(181.6*(1+(0.03*iden)+0.062*(t/303)*(iden**2.5)))/((lambd-0.634-3*iden)*exp(4.18*((t-303)/t)))
        rm_model=(pm_model*l)/a
        nohm_model=i*(rm_model+rc)
        vcell_model=enernst-nact_model-nohm_model-nconc_model
        #print(enernst)
        #print(nact_model)
        #print(nohm_model)
        #print(nconc_model)
        vstack_model=ns*vcell_model
        #print(vstack_model)             
        nact_actual=-(-0.944957+(0.00301801*t)+(0.00007401*t*log1p(co2))+(-0.000188*t*log1p(i)))
        #print(-0.944957,0.00301801*t,0.00007401*t*log1p(co2),-0.000188*t*log1p(i))
        nconc_actual=-0.02914489*log1p(1-(iden/ilimitden))
        pm_actual=(181.6*(1+(0.03*iden)+0.062*(t/303)*(iden**2.5)))/((23-0.634-3*iden)*exp(4.18*((t-303)/t)))
        rm_actual=(pm_model*l)/a
        nohm_actual=i*(rm_model+0.0001)
        vcell_actual=enernst-nact_actual-nohm_actual-nconc_actual
        vstack_actual=ns*vcell_actual
        #print(enernst,nact_actual,nohm_actual,nconc_actual)
        #print(enernst,nact_model,nohm_model,nconc_model)
        #print(vstack_actual)
        sum+=(vstack_actual-vstack_model)**2
    #print(sum)
    return sum


def pop(itr=50):
    lst=[]
    
    for x in range(itr):
        innerlst=[]
        x1=individual_x1(-0.952,-0.944)
        x2=individual_x2(0.003,0.003)
        x3=individual_x3(0.000074,0.000078)
        x4=individual_x4(-0.000198,-0.000188)
        lambd=individual_lambda(14,23)
        rc=individual_rc(0.0001,0.0008)
        b=individual_b(0.055,0.055)
        innerlst.extend([x1,x2,x3,x4,lambd,rc,b])
        lst.append(innerlst)
    #print(lst)
    return lst

def evolve(retain=0.6,random_select=0.05,mutate=0.01):
    #selection
    lst=pop()
    graded=[x for x in sorted(lst)]
    retain_length=int(len(graded)*retain)
    parents=graded[:retain_length]
    for individual in graded[retain_length:]:
        if(random_select)>random():
            parents.append(individual)

    #crossover parents to create children
    parents_length=len(parents)
    #print(parents_length)
    desired_length=len(lst)-parents_length
    #print(desired_length)
    children=[]
    while len(children)<desired_length:
        male=randint(0,parents_length-1)
        female=randint(0,parents_length-1)
        if male!=female:
            male=parents[male]
            #print(male)
            female=parents[female]
            #print(female)
            half=int(len(male)/2)
            #print(half)
            child=male[:half]+female[half:]
            children.append(child)
    parents.extend(children)

    #mutate some individuals
    for individual in parents:
        if mutate>random():
            pos_to_mutate=randint(0,len(individual)-1)
            individual[pos_to_mutate]+=0.01
    #print(parents)        
    return parents

def run_avg_fitness():
    parents=evolve()
    #print(parents)
    sum_fitness=0
    for x in range(len(parents)):
        x1=parents[x][0]
        #print(parents[x][0])
        x2=parents[x][1]
        #print(parents[x][1])
        x3=parents[x][2]
        #print(parents[x][2])
        x4=parents[x][3]
        #print(parents[x][3])
        lambd=parents[x][4]
        #print(parents[x][4])
        rc=parents[x][5]
        #print(parents[x][5])
        b=parents[x][6]
        #print(parents[x][6])
        #print(fitness(x1,x2,x3,x4,lambd,rc,b))
        sum_fitness+=fitness(x1,x2,x3,x4,lambd,rc,b)
    #print(sum_fitness)
    #print(len(parents))
    return (sum_fitness/len(parents))

def main():
    for i in range(-15,15,1):
        seed(a=i)
        best_run=run_avg_fitness()
        avg=0
        for x in range(10000):
            avg=run_avg_fitness()
            if avg<best_run:
                best_run=avg
        print(best_run)

main()

