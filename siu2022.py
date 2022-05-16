# This code is written by Rahmi Elibol at February 2022 for siu2022.


from brian2 import *
start_scope()


# Parameters
#neuron parameters
Vthr = 30 * mvolt
EL = -75 * mV


#synaptic parameters
tau_s = 1 * ms
we = 0.1 *amp/mV
wi = 0.1 *amp/mV
Vi = -90 * mV
Ve = 0 * mV
dly=(3+rand())*ms


### synaptic weigths
w_cse = 1.00 * we        

w_vse = 2.00 *we
w_vsi = 1.00 *wi



par_percent=10


#number_of_neurons_in_cortex=900
number_of_neurons_in_P = 100
number_of_neurons_in_S1 = 2
number_of_neurons_in_S2 = 2
number_of_neurons_in_D = 100




#############   ---- Poisson Groups ------ #####################

PG_P_D_1011 = PoissonGroup(number_of_neurons_in_P, 10 * Hz)  #nominal value: 5*Hz


print('Equations')


eqs_dyn = """
dv/dt=(0.04/ms/mV)*v**2+(5/ms)*v+140*mV/ms-u/ms+I*mV/(amp*ms)+Is*mV/(amp*ms) : volt
du/dt=a*(b*v-u)/ms                                           : volt
I : amp
Is=ge*(Ve-v)+gi*(Vi-v) :  amp
dge/dt=-ge/tau_e	: amp/volt
dgi/dt=-gi/tau_i	: amp/volt
a : 1
b : 1
c : volt
d : volt
tau_e : second
tau_i : second
"""


tau_glu=2*ms
tau_DA=1.5*ms
tau_i_msn=tau_s

eqs_s = """
dv/dt=(0.04/ms/mV)*v**2+(5/ms)*v+140*mV/ms-u/ms+I*mV/(amp*ms)+Is*mV/(amp*ms) : volt
du/dt=a*(b*v+k*mV-u)/ms                                           : volt
I : amp
I_Glu=g_glu*(Ve-v) :  amp
I_DA=g_DA*(V_DA-v) :  amp
I_Ach=g_Ach*(Vi-v)  :  amp
I_GABA=g_GABA*(Vi-v) :  amp
Is=I_Glu+I_DA+I_Ach+I_GABA :  amp
dg_glu/dt=-g_glu/tau_glu	: amp/volt
dg_DA/dt=-g_DA/tau_DA	: amp/volt
dg_Ach/dt=-g_Ach/tau_i_msn	: amp/volt
dg_GABA/dt=-g_GABA/tau_i_msn	: amp/volt
a : 1
b : 1
c : volt
d : volt
k : 1
V_DA : volt
"""




eqs_reset = '''
v = c
u = u+d
'''


###################################################
####### ----------- Groups -----#######
###################################################


P = NeuronGroup(number_of_neurons_in_P, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)


for i in range(number_of_neurons_in_P):
    P.a[i] = 0.02*((100-par_percent+2*par_percent*rand())/100)
    P.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    P.c[i] = -65*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    P.d[i] = 8*(100-par_percent+2*par_percent*rand())/100* mvolt 
    P.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    P.u[i] = (-14.5*((100-par_percent+2*par_percent*rand())/100))*mvolt
P.tau_e = tau_s
P.tau_i = tau_s
    
    



S1= NeuronGroup(number_of_neurons_in_S1, model=eqs_s, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_S1):
    S1.a[i] = 0.02*((100-par_percent+2*par_percent*rand())/100)
    S1.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    S1.c[i] = -56*((100-par_percent+2*par_percent*rand())/100) * mvolt      # NV:-52  
    S1.d[i] = 0.4*((100-par_percent+2*par_percent*rand())/100) * mvolt     #NV: 1.9
    S1.v[i] = (EL-15*mV)*((100-par_percent+2*par_percent*rand())/100)
    S1.u[i] = 35*((100-par_percent+2*par_percent*rand())/100)*mV
    S1.k[i] = 35*((100-2*par_percent+4*par_percent*rand())/100)
S1.V_DA = 0*mV    
    


S2= NeuronGroup(number_of_neurons_in_S2, model=eqs_s, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_S2):
    S2.a[i] = 0.02*((100-par_percent+2*par_percent*rand())/100)
    S2.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    S2.c[i] = -55*((100-par_percent+2*par_percent*rand())/100) * mvolt        #NV: -52
    S2.d[i] = 0.4*((100-par_percent+2*par_percent*rand())/100) * mvolt       #NV: 1.9   
    S2.v[i] = (EL-15*mV)*((100-par_percent+2*par_percent*rand())/100)
    S2.u[i] = 25*((100-par_percent+2*par_percent*rand())/100)*mV
    S2.k[i] = 20*((100-2*par_percent+4*par_percent*rand())/100)
S2.V_DA = -90*mV



D = NeuronGroup(number_of_neurons_in_D, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)


for i in range(number_of_neurons_in_D):
    D.a[i] = 0.02*((100-par_percent+2*par_percent*rand())/100)
    D.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    D.c[i] = -70*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    D.d[i] = 8*((100-par_percent+2*par_percent*rand())/100) * mvolt 
    D.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    D.u[i] = -14.5 *((100-par_percent+2*par_percent*rand())/100) *mV
D.tau_e = tau_s
D.tau_i = tau_s




    


P.I = 0*amp
S1.I = 0*amp
S2.I = 0*amp
D.I = 0*amp

#==============================================================================
#==============================================================================
#############################
print('Synapses')
#############################
#############################

S01_1011_111 = Synapses(PG_P_D_1011, P,  'w :siemens', delay=dly, on_pre='ge += w')
S01_1011_111.connect(True, p = 0.25)

# # STDP
# stdp parameters
#==============================================================================
 
taupre = 20*ms
taupost = taupre
 
w_base=80
gmax = 150
dApre = .01
dApost = -dApre * 1.05
dApost *= gmax
dApre *= gmax


S21_111_123 = Synapses(P, S1, 
              '''w : 1
                 dApre/dt = -Apre / taupre : 1 (event-driven)
                 dApost/dt = -Apost / taupost : 1 (event-driven)''',
              on_pre='''g_glu += w*siemens
                     Apre += dApre
                     w = clip(w + Apost, 0, gmax)''',
              on_post='''Apost += dApost
                      w = clip(w + Apre, 0, gmax)''',
              )
S21_111_123.connect(True)
 
 #initial value for weights
 
 
S21_111_123.w = 'rand() * 20 + w_base'

print("Number of synapses from P: "+str(S21_111_123.N)+" average:"+str(S21_111_123.N/number_of_neurons_in_S1))


S24_131_123 = Synapses(D, S1,  'w :siemens', delay=dly, on_pre='g_DA += w')
S24_131_123.connect(True)
S24_131_123.w=w_vse


print("Number of synapses from D : "+str(S24_131_123.N)+" average:"+str(S24_131_123.N/number_of_neurons_in_S1))
print("----------------------------------------------")



S29_111_124 = Synapses(P, S2, 
              '''w : 1
                 dApre/dt = -Apre / taupre : 1 (event-driven)
                 dApost/dt = -Apost / taupost : 1 (event-driven)''',
              on_pre='''g_glu += w*siemens
                     Apre += dApre
                     w = clip(w + Apost, 0, gmax)''',
              on_post='''Apost += dApost
                      w = clip(w + Apre, 0, gmax)''',
              )
S29_111_124.connect(True)
 
#initial value for weights
 
S29_111_124.w = 'rand() * 20 + w_base'

print("Number of synapses from P : "+str(S29_111_124.N)+" average:"+str(S29_111_124.N/number_of_neurons_in_S1))


S32_131_124 = Synapses(D, S2,  'w :siemens', delay=dly, on_pre='g_DA += w')
S32_131_124.connect(True)
S32_131_124.w=w_vsi

print("Number of synapses from D : "+str(S32_131_124.N)+" average:"+str(S32_131_124.N/number_of_neurons_in_S1))


S54_1011_131 = Synapses(PG_P_D_1011, D,  'w :siemens', delay=dly, on_pre='ge += w')
S54_1011_131.connect(True, p = 0.25)
S54_1011_131.w=we*0.5

#==============================================================================
#==============================================================================




import time
init_time=time.time()

################################################################
####################------------ First 100ms --------- #########
################################################################

duration1=100*ms


S01_1011_111.w=we*0
S54_1011_131.w=we*0


print("----------------------------------------------")
print('100 ms initial conditions delay')
print("sim_time="+str(duration1))
run(duration1, report='text')




S01_1011_111.w=we    ############# ACA Cortex
S54_1011_131.w=we     ############# VTA DA



#==============================================================================
#==============================================================================

#==============================================================================
#==============================================================================


#==============================================================================
#==============================================================================


###############################################################################
##################### ---- Monitors ---- ######################################
###############################################################################
print('Monitors')




trace_P = StateMonitor(P, 'v', record=0)
spikes_P = SpikeMonitor(P)


trace_S1 = StateMonitor(S1, 'v', record=True)
spikes_S1 = SpikeMonitor(S1)
trace_I_s_msnd1_shell = StateMonitor(S1, 'Is', record=True)


trace_S2 = StateMonitor(S2, 'v', record=True)
spikes_S2 = SpikeMonitor(S2)
trace_I_s_msnd2_shell = StateMonitor(S2, 'Is', record=True)

trace_D = StateMonitor(D, 'v', record=0)
spikes_D = SpikeMonitor(D)

###############################################################################
###############################################################################
## --------- Synaptic Monitors ----------------

mon_S21_111_123 = StateMonitor(S21_111_123, 'w', record=True)

mon_S29_111_124 = StateMonitor(S29_111_124, 'w', record=True)

###############################################################################
###############################################################################

###############################################################################
print("----------------------------------------------")
print("start")

number_of_scenario=2


#########################  ----- Scenarios   --------------  ####################

##-------------------- Scenario 1 --------------------------------##
#### for testing the code


if number_of_scenario==1:

    print("----------------------------------------------")
    print('Scenario 1')
    run(100*ms,report='text')



elif number_of_scenario==2:
    
    print("----------------------------------------------")
    print('Scenario 2')
        
    print("Scenario 2: 1 : P var")
    
    S01_1011_111.w=we    ############# ACA Cortex
    S54_1011_131.w=we*0.3     ############# VTA DA
    
    duration=1000*ms
    
    run(duration,report='text')
        
        
    print("P :"+str(spikes_P.num_spikes)+"  ---->   : " +str(spikes_P.num_spikes/number_of_neurons_in_P))
    print("S1 : "+str(spikes_S1.num_spikes)+"  ---->  : " +str(spikes_S1.num_spikes/number_of_neurons_in_S1))
    print("S2 : "+str(spikes_S2.num_spikes)+"  ---->  : " +str(spikes_S2.num_spikes/number_of_neurons_in_S2))
    print("D             : "+str(spikes_D.num_spikes)+"  ---->       : " +str(spikes_D.num_spikes/number_of_neurons_in_D))

    print("----------------------------------------------")
    print("Scenario 2: 2 : D var")
    
    S01_1011_111.w=we*0.3    ############# ACA Cortex
    S54_1011_131.w=we     ############# VTA DA
    

    run(duration,report='text')

        
        
    print("P :"+str(spikes_P.num_spikes)+"  ---->   : " +str(spikes_P.num_spikes/number_of_neurons_in_P))
    print("S1 : "+str(spikes_S1.num_spikes)+"  ---->  : " +str(spikes_S1.num_spikes/number_of_neurons_in_S1))
    print("S2 : "+str(spikes_S2.num_spikes)+"  ---->  : " +str(spikes_S2.num_spikes/number_of_neurons_in_S2))
    print("D             : "+str(spikes_D.num_spikes)+"  ---->       : " +str(spikes_D.num_spikes/number_of_neurons_in_D))

    
    print("----------------------------------------------")
    print("Scenario 2: 3 : P ve D var")
    
    S01_1011_111.w=we    ############# ACA Cortex
    S54_1011_131.w=we     ############# VTA DA
    
    
    run(duration,report='text')

        
        
    print("P :"+str(spikes_P.num_spikes)+"  ---->   : " +str(spikes_P.num_spikes/number_of_neurons_in_P))
    print("S1 : "+str(spikes_S1.num_spikes)+"  ---->  : " +str(spikes_S1.num_spikes/number_of_neurons_in_S1))
    print("S2 : "+str(spikes_S2.num_spikes)+"  ---->  : " +str(spikes_S2.num_spikes/number_of_neurons_in_S2))
    print("D             : "+str(spikes_D.num_spikes)+"  ---->       : " +str(spikes_D.num_spikes/number_of_neurons_in_D))

    

else:
    print('No scenario !!!!')


    
print("----------------------------------------------")
filename=str(time.time())
final_time=time.time()
print("Total simulation time:",str(final_time-init_time))

#### ------------------------------------------------------------------- ####
#### ------------------------------------------------------------------- ####
#### ------------------------------------------------------------------- ####


print("P : "+str(spikes_P.num_spikes)+"  ---->  : " +str(spikes_P.num_spikes/number_of_neurons_in_P))
print("S1     : "+str(spikes_S1.num_spikes)+"  ---->  : " +str(spikes_S1.num_spikes/number_of_neurons_in_S1))
print("S2     : "+str(spikes_S2.num_spikes)+"  ---->  : " +str(spikes_S2.num_spikes/number_of_neurons_in_S2))
print("D             : "+str(spikes_D.num_spikes)+"  ---->  : " +str(spikes_D.num_spikes/number_of_neurons_in_D))



print("----------------------------------------------")





#########################################
############# ----- Weights --------
#########################################
print("------------- weights  ------")

dec_inc_S1=[]
max_temp=0
min_temp=0

for i in range(len(mon_S21_111_123.w)):
    array_temp=[]
    max_temp=0
    min_temp=0
    array_temp=mon_S21_111_123.w[i]
    max_temp=max(array_temp)
    min_temp=min(array_temp)
    max_index_k=where(array_temp==max_temp)
    max_index=max_index_k[0][0]   
    min_index_k=where(array_temp==min_temp)
    min_index=min_index_k[0][0]
    if max_index<min_index:
        dec_inc_S1.append(min_temp-max_temp)
    elif max_index>min_index:
        dec_inc_S1.append(max_temp-min_temp)
    else:
        dec_inc_S1.append(max_temp-min_temp)


max_index_S1=dec_inc_S1.index(max(dec_inc_S1))
min_index_S1=dec_inc_S1.index(min(dec_inc_S1))


dec_inc_S2=[]
max_temp=0
min_temp=0

for i in range(len(mon_S29_111_124.w)):
    array_temp=[]
    max_temp=0
    min_temp=0
    array_temp=mon_S29_111_124.w[i]
    max_temp=max(array_temp)
    min_temp=min(array_temp)
    max_index_k=where(array_temp==max_temp)
    max_index=max_index_k[0][0]   
    min_index_k=where(array_temp==min_temp)
    min_index=min_index_k[0][0]
    if max_index<min_index:
        dec_inc_S2.append(min_temp-max_temp)
    elif max_index>min_index:
        dec_inc_S2.append(max_temp-min_temp)
    else:
        dec_inc_S2.append(max_temp-min_temp)


max_index_S2=dec_inc_S2.index(max(dec_inc_S2))
min_index_S2=dec_inc_S2.index(min(dec_inc_S2))



figure()
subplot(211)
plot(mon_S21_111_123.t/ms,mon_S21_111_123.w[max_index_S1],'k')
plot(mon_S21_111_123.t/ms,mon_S21_111_123.w[min_index_S1],'k')
plot(mon_S29_111_124.t/ms,mon_S29_111_124.w[max_index_S2],'r')
plot(mon_S29_111_124.t/ms,mon_S29_111_124.w[min_index_S2],'r')


weights_change_S1=[]
weights_change_S2=[]

if number_of_scenario==2:
    total_step=60
    step_time=len(mon_S21_111_123.t)/total_step
else:
    total_step=len(mon_S21_111_123.t)
    step_time=len(mon_S21_111_123.t)/len(mon_S21_111_123.t)


print("total_steps_synaptic_weights : "+str(total_step))
print("step_time_synaptic_weights : "+str(step_time))

for t in range(total_step):
    print(int(t*step_time))
    print(100*step_time*t/len(mon_S21_111_123.t))
    total_temp_4_t_S1=0
    total_temp_4_t_S2=0
    for i in range(len(mon_S21_111_123.w)):
        total_temp_4_t_S1=total_temp_4_t_S1+mon_S21_111_123.w[i][int(t*step_time)]
        total_temp_4_t_S2=total_temp_4_t_S2+mon_S29_111_124.w[i][int(t*step_time)]
    weights_change_S1.append(total_temp_4_t_S1)
    weights_change_S2.append(total_temp_4_t_S2)


a_S1=weights_change_S1[-1]/weights_change_S1[0]
print("area s1 ="+str(a_S1))
a_S2=weights_change_S2[-1]/weights_change_S2[0]
print("area s2 ="+str(a_S2))

y_S1=((weights_change_S1[-1]/S21_111_123.N)-(weights_change_S1[0]/S21_111_123.N))
print("y s1 ="+str(y_S1))
y_S2=((weights_change_S2[-1]/S29_111_124.N)-(weights_change_S2[0]/S29_111_124.N))
print("y s2 ="+str(y_S2))


y5=weights_change_S1
y6=weights_change_S2


print("wS1 total increase in delta t1 : "+ str(100*(y5[20]-y5[0])/y5[0]))
print("wS1 total increase in delta t2: "+ str(100*(y5[40]-y5[20])/y5[20]))
print("wS1 total increase in delta t3: "+ str(100*(y5[-1]-y5[40])/y5[40]))
print("wS1 total increase : "+ str(100*(y5[-1]-y5[0])/y5[0]))

print("wS2 total increase in delta t1: "+ str(100*(y6[20]-y6[0])/y6[0]))
print("wS2 total increase in delta t2: "+ str(100*(y6[40]-y6[20])/y6[20]))
print("wS2 total increase in delta t3: "+ str(100*(y6[-1]-y6[40])/y6[40]))
print("wS2 total increase : "+ str(100*(y6[-1]-y6[0])/y6[0]))


subplot(212)
plot(weights_change_S1)
plot(weights_change_S2)



#########################################
############# --- Figures --- #########
#########################################

figure()
subplot(411)
plot(trace_P.t / ms, trace_P[0].v / mV)
ylabel('P')


subplot(412)
plot(trace_D.t / ms, trace_D[0].v / mV)
ylabel('D')


subplot(413)
plot(trace_S1.t / ms, trace_S1[0].v / mV)
ylabel('S1')



subplot(414)
plot(trace_S2.t / ms, trace_S2[0].v / mV)
xlabel('time, ms')
ylabel('S2')


show()

print("----------------------------------------------")
print(" - - -     The End     - - - ")
print("----------------------------------------------")
