function [p_received,db_received]=received_pressure(p_vector,TL_matrix)

% function to calculate the received pressure at a set of recorders
% by Sebastian Menze, 2015
%
% inputs: source pressure vector and the transmission loss matrix connecting
% sources and receivers
%
% outpus: received pressure in uPa dB and received level in dB re 1 uPa

for i_receiver=1:size(TL_matrix,2)
       
    for i_source=1:size(TL_matrix,1)
       
        single_source_RL= 20*log10(p_vector(i_source)) + TL_matrix(i_source,i_receiver) ;
        if single_source_RL < 0
             received_db(i_source)=0 ;
        else
             received_db(i_source)=single_source_RL ;
        end     
            clear single_source_RL  
    end  
    
      p_received(i_receiver)=sum(10.^(received_db/20));
      db_received(i_receiver)=  20*log10(p_received(i_receiver));
     clear received_db
end


end