# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 17:32:49 2018

@author: Seb
"""

import os
import glob
import multiprocessing  
import scipy.io
import numpy as np
import random
import time
from functools import partial

def received_pressure(p_vector,tl_db):
    p_received=list()
    for i_receiver in range(0,tl_db.shape[1]):  
        received_db=list()
        for i_source in range(0,tl_db.shape[0]):       
            single_source_RL= 20*np.log10(p_vector[i_source]) + tl_db[i_source,i_receiver] 
            if single_source_RL < 0:
                 received_db.append(0) 
            else:
                received_db.append(single_source_RL)        
        p_received.append( sum(10**(np.array(received_db)/20)) )
    db_received =20*np.log10(np.array(p_received))      
    return db_received


def parafunc(i_solution,sim_whale_location_old,p_sim_whale,tl_db,db_true,n_sim_whale,n_source_locations,sigma_db,n_iterations):
    
    likelihood=list()
    sse=list()
    t_start=1;
    t_exponent=2
    tic = time.clock()
    p_sim_whale_par=p_sim_whale[i_solution]
    print('Solution '+str(i_solution)+' running!')
    # serial loop through iterations
    for i_iteration in range(n_iterations):
        
        
        # update temperature
        temp= t_start*np.power( (1-(i_iteration)/n_iterations) , t_exponent )
        
        ix_selected_simwhale=random.randint(0,n_sim_whale-1)
        new_location=random.randint(1,n_source_locations)
        sim_whale_location_new=sim_whale_location_old
        sim_whale_location_new[ix_selected_simwhale]=new_location
        
        if not 'likelihood_old' in locals():
            pressure_old=np.empty(n_source_locations)
            for i_hexagonal_source in range(1,n_source_locations):
                pressure_old[i_hexagonal_source]=sum(sim_whale_location_old==i_hexagonal_source)*p_sim_whale_par
            db_old=received_pressure(pressure_old,tl_db)
            # old misfit and likelihood
            sse_old = 0.5 * np.sum(  ((db_old - db_true)**2) / (sigma_db**2)   )
            likelihood_old= np.exp( - sse_old ) 
        
        pressure_new=np.empty(n_source_locations)
        for i_hexagonal_source in range(1,n_source_locations):
            pressure_new[i_hexagonal_source]=sum(sim_whale_location_new==i_hexagonal_source)*p_sim_whale_par 
        db_new=received_pressure(pressure_new,tl_db)
        
        # new misfit and likelihood
        sse_new = 0.5 * np.sum(  ((db_new - db_true)**2) / (sigma_db**2)   )
        likelihood_new= np.exp( - sse_new ) 
        
        if likelihood_new==0:    
            if sse_new <=  sse_old:
                sim_whale_location_old[ix_selected_simwhale]=new_location  
                likelihood_old=likelihood_new
                sse_old=sse_new
            
            else:    
                if  random.expovariate(1/temp) >1:            
                    sim_whale_location_old[ix_selected_simwhale]=new_location  
                    likelihood_old=likelihood_new
                    sse_old=sse_new
        else:
             if likelihood_old <= likelihood_new:
                sim_whale_location_old[ix_selected_simwhale]=new_location  
                likelihood_old=likelihood_new
                sse_old=sse_new
             else:    
                if  random.expovariate(1/temp) >1:            
                    sim_whale_location_old[ix_selected_simwhale]=new_location        
                    likelihood_old=likelihood_new
                    sse_old=sse_new
        
        likelihood.append(likelihood_old)
        sse.append(sse_old)

#        print('Iteration '+str(i_iteration).zfill(6)+' L(m): '+str(likelihood_old)+' time: '+str(toc-tic) )
    toc = time.clock()    
    print('solution '+str(i_solution)+' completed in: '+str(toc-tic)+' L(m)= '+str(likelihood_old))    
    return sim_whale_location_old,likelihood,sse


if __name__ == '__main__':


#    workfolder = r'C:\Users\Seb\Documents\passive_acoustic_work\weddell_sea_scenarios'
#    os.chdir(workfolder)
    n_iterations=20000
    n_solutions=50
    p_min=1e11
    p_max=1e13
    t_exponent=2
    sigma_db=1
    
    scenariofiles=glob.glob('scenario*.mat')
    random.shuffle(scenariofiles)
    
    for selectedfile in scenariofiles:

        if not os.path.isfile('solution'+selectedfile[8:-3]+'mat'):
            print('--> File: '+selectedfile)
            m=scipy.io.loadmat(selectedfile)
            
            n_sim_whale=np.max( m['sim_sources']['id'][0,0] )
            n_source_locations=np.max( m['sim_sources']['id'][0,0] )
            tl_db=m['tl_db']
            db_true=m['recorder']['db_received'][0,0]
        #    t1 = time.clock()
        
            startlocations=np.array([random.randint(1,n_source_locations) for i in range(n_source_locations)])
            p_sim_whale=np.linspace(p_min,p_max,num=n_solutions)
            p_sim_whale=p_sim_whale/n_sim_whale;
            
            cpucounts=multiprocessing.cpu_count()
            
            print('Working with '+str(cpucounts)+' cores')
            pool = multiprocessing.Pool(processes=cpucounts)
            
            i_solutions=range(n_solutions)
            a=pool.map(partial(parafunc,sim_whale_location_old=startlocations,p_sim_whale=p_sim_whale,tl_db=tl_db,db_true=db_true,n_sim_whale=n_sim_whale,n_source_locations=n_source_locations,sigma_db=sigma_db,n_iterations=n_iterations), i_solutions)
            pool.close      
        #    t2 = time.clock()
                
            likelihood=np.empty(n_solutions)
            sse=np.empty(n_solutions)
            est_p_sum=np.empty(n_solutions)
            est_p=np.empty([n_solutions,n_sim_whale])
            likelihood_iterations=np.empty([n_solutions,n_iterations])
            sse_iterations=np.empty([n_solutions,n_iterations])
            
            for i in range(n_solutions):
                for j in range(n_sim_whale):
                     #est_p[i,j]=sum( np.in1d(a[i][0] , j) ) * p_sim_whale[i]
                     est_p[i,j]=sum( np.array(a[i][0])==np.array(j) ) * p_sim_whale[i]
                     likelihood_iterations[i,:]=a[i][1]
                     sse_iterations[i,:]=a[i][2]
                     likelihood[i]=a[i][1][-1]
                     sse[i]=a[i][2][-1]
                     est_p_sum[i]=p_sim_whale[i]*n_sim_whale
        
            scipy.io.savemat('solution'+selectedfile[8:-3]+'mat', mdict={'likelihood': likelihood,'sse':sse,'est_p_sum':est_p_sum,'est_p':est_p,'likelihood_iterations':likelihood_iterations,'sse_iterations':sse_iterations})


    
    