function sa_parameter_estimation(n_iterations,n_solutions,t_exponent,p_min,p_max,sim_source_location,recorder_location,true_rl_file,TL_db,sigma_db,target_file)

% SIMULATED ANNEALING PARAMETER ESTIMATION FUNCETION
% by Sebastian Menze, 2015
% used to invert fin whale chorus received level observations 
% to most likely source grid pattern
% 
% inputs: 
%   number of iterations, number of solutions, cooling exponent,
%   minimum total pressure, maxiumim total pressure, pathto soiulated source
%   grid file, path to recorder location file, path to received level file,
%   transmission loss matrix, forward model standard deviation, name of
%   outpit fie
%
% outputs:
%   solutions matrix in uPa [columns: number of solution , rows: number of source grid nodes]  
%   solutions matrix in dB re 1 uPa [columns: number of solution , rows: number of source grid nodes]  
%   vector from 1 : number of iterations
%   misfit matrix [columns: number of solutions, rows: number of iterations]
%   likelihood matrix [columns: number of solutions, rows: number of iterations]

% read source grid node locations
[sim_source_id,sim_source_lat,sim_source_lon] = import_loc_file(sim_source_location) ;
n_sim_whale=sim_source_id(end);
n_source_locations=sim_source_id(end);

% read recorder positions
% recorder_location='recorder_location_files/recorders_geodesic_47.loc';
[recorders_id,recorders_lat,recorders_lon] = import_loc_file(recorder_location) ;

% read true rl 
% true_rl_file='true_source_files/tl_cylindrical_ts_gauss_r47.rl';
[true_rl_id,true_rl_p,true_rl_db] = import_true_rl(true_rl_file) ;

% place one simulated source on each node and generate homogenous a priori source
% grids between the minimum and maximum total pressure
sim_whale_location_old_all=repmat(sim_source_id',n_solutions,1) ;
if p_min<0
    p_min=0;
end
p_sim_whale=linspace(p_min,p_max,n_solutions);
p_sim_whale=p_sim_whale/n_sim_whale;

%%% use simulated annealing
t_start=1;

for i_iteration=1:n_iterations
    
    tic
    % update temperature
    temp= t_start*(1-(i_iteration)./n_iterations).^t_exponent;
    mu=temp;
    
    sim_whale_location_old_all_new=[];
    p_sim_whale_par=[];
        p_sim_whale_par_all=[];

    parfor i_solution=1:n_solutions
    
 
  %%% move random sim_whale to random new loation
  ix_selected_sim_whale=randi(n_sim_whale);
  new_location=randi(n_source_locations);
  
   sim_whale_location_old=sim_whale_location_old_all(i_solution,1:n_sim_whale);

%  clear sim_whale_location_new
  sim_whale_location_new=sim_whale_location_old;
  sim_whale_location_new(ix_selected_sim_whale)=new_location;
  
  %%% calculate new pressure
p_sim_whale_par=p_sim_whale(i_solution);
[pressure_new]=parfor_pressure_helper(n_source_locations,sim_whale_location_new,p_sim_whale_par)
[pressure_old]=parfor_pressure_helper(n_source_locations,sim_whale_location_old,p_sim_whale_par)

% new misfit and likelihood
[p_received,db_received]=received_pressure(pressure_new,TL_db) ;
 sse_new = 0.5 * sum(  ((db_received' - true_rl_db).^2) ./ (sigma_db^2)   );  
likelihood_new= exp( - sse_new ) ;

 % old misfit and likelihood
[p_received,db_received]=received_pressure(pressure_old,TL_db) 
 sse_old = 0.5 * sum(  ((db_received' - true_rl_db).^2) ./ (sigma_db^2)   );  
likelihood_old= exp( - sse_old ) ;

if likelihood_new==0
   if sse_new <=  sse_old
% move sim_whale
       sim_whale_location_old(ix_selected_sim_whale)=new_location;        
   else    
       if  random('exp',mu)>1            
       sim_whale_location_old(ix_selected_sim_whale)=new_location;
       end        
   end             
else
   if likelihood_old <= likelihood_new
% move sim_whale
       sim_whale_location_old(ix_selected_sim_whale)=new_location;        
   else    
       if  random('exp',mu)>1            
       sim_whale_location_old(ix_selected_sim_whale)=new_location;
       end        
   end
end   
     
     % update simulates source locations and likelihood values due to
     % parfor constraints
     performance_sse(i_iteration,i_solution)=sse_old;      
     performance_likelihood(i_iteration,i_solution)=likelihood_old;  
     sim_whale_location_old_all_new= [sim_whale_location_old_all_new; sim_whale_location_old ];
     p_sim_whale_par_all=[p_sim_whale_par_all,p_sim_whale_par];
      
    end
    
    p_sim_whale=p_sim_whale_par_all;
    sim_whale_location_old_all=sim_whale_location_old_all_new;
    performance_iteration(i_iteration)=i_iteration;


disp(['done: ',num2str(i_iteration/n_iterations*100),'% and best L(m): ',num2str((max(performance_likelihood(i_iteration,:)))),' min misfit: ',num2str((min(performance_sse(i_iteration,:)))),' dB'])
toc
    
end

% calculate source grid based on the source locations obtained from
% simulated annealing
    for i_solution=1:n_solutions
for i_hexagonal_source=1:n_source_locations
solutions_p(i_solution,i_hexagonal_source)=sum(sim_whale_location_old_all(i_solution,:)==i_hexagonal_source)*p_sim_whale(i_solution) ;
end  
end

 solutions_db=20.*log10((solutions_p));
 solutions_db(solutions_db<0)=0;
 
 % save results as mat file
 save(target_file,'solutions_p','solutions_db','performance_iteration','performance_likelihood','performance_sse','p_sim_whale')

end
    