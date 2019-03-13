function [tl_db]=tl_matrix_cylindrical_spreading2(sim_sources,recorders,f_khz,waterdepth)

% function to calculate TL matrix using the cylindrical spreading law 
% by Sebastian Menze, 2015
%
% input: 
% path to simulated source grid file, path to recorder location
% file, file name under wich TL matrix will be stored, frequency in kHz,
% range until spherical spreading is assumed (average water depth)
%
% output:
% TL matrix in dB re 1 upa as .mat file
%
% example:
% tl_cylindrical_spreading('sim_sources.location_files/sim_sources_geodesic_195.loc','recorder_location_files/recorders.geodesic_47.loc','ss195_r47.mat',0.2,100)

% waterdepth=4000
% f=0.2;
absorption_factor = sound_abs(2,35,10,f_khz,7.5);

% absorption_factor=10^-3
% db_uncertainty=2
nsize=numel(sim_sources.lat);

l=round(linspace(1,nsize,100));

for i_source=1:nsize
 
if  ismember(i_source,l)
    disp( [ num2str(i_source/nsize*100),' percent done!'] ) 
end

for i_receiver=1:numel(recorders.lat)
    
distance_in_km(i_source,i_receiver)=deg2km(distance(sim_sources.lat(i_source),sim_sources.lon(i_source),recorders.lat(i_receiver),recorders.lon(i_receiver))) ;

if distance_in_km(i_source,i_receiver) > 0
    
    if distance_in_km(i_source,i_receiver)*1000 < waterdepth
        
   % spherical loss until water depth  
TL_spreading=20*log10( distance_in_km(i_source,i_receiver) *1000 ) ;
tl_db(i_source,i_receiver)= absorption_factor * distance_in_km(i_source,i_receiver) + TL_spreading ;

    else
  % spherical nd cylindrical loss u
TL_spreading= 20*log10(waterdepth) +  10*log10( distance_in_km(i_source,i_receiver)*1000 /waterdepth  ) ;
tl_db(i_source,i_receiver)= absorption_factor * distance_in_km(i_source,i_receiver) + TL_spreading ;
       
    end
   
else
 tl_db(i_source,i_receiver)= 0;   

end
clear TL_spreading

end                
end
tl_db=-tl_db;
% save(target_filename,'tl_db')

