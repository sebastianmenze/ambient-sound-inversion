%% Supplementary code for the paper "Estimating the spatial distribution of 
% vocalizing animals from ambient sound spectra using widely spaced
% recorder arrays and inverse modelling"

% by Sebastian Menze, Daniel Zitterbart, Martin Biuw and Olaf Boebel
% submitted to JASA in 2019
% contact: sebastian.menze@hi.no
% created and tested with Matlab 2016a

%% these liberaries needs to be downloaded separately:

% Grid sphere package to generate geodesic grids:
% https://se.mathworks.com/matlabcentral/fileexchange/28842-grid-sphere
% 
% cmocean colormaps:
% https://se.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
% 
% m_map mapping package:
% https://www.eoas.ubc.ca/~rich/map.html
% 
% freezeColors code to use multiple colormaps on the same figure axis:
% https://se.mathworks.com/matlabcentral/fileexchange/7943-freezecolors-unfreezecolors
% 
% function to calculate seawater sound absoption such as this one:
% https://se.mathworks.com/matlabcentral/fileexchange/28653-acoustic-absorption-in-seawater
% 
% function to generate random matrices:
% https://se.mathworks.com/matlabcentral/fileexchange/59305-randnd

%% Contents of this folder:

% this matlab script to generate and analyse inversion test scenarios
% codeforpaper.m
% 
% parallel parameter optimization pyhton code (simulated annealing):
% paraanneal.py
% 
% shell script to compute MCMCs on unix cluster:
% paraanneal.sh
% 
% matlab function to calculate geometrical transmission loss:
% tl_matrix_geometrical_spreading.m
% 
% matlab function to calculate RL:
% received_pressure.m
% 
% matlab code for parameter optimization (in older format, not used for paper)
% sa_parameter_estimation.m

%% add liberaries

addpath(genpath('m_map1.4'))
addpath(genpath('GridSphere'))
addpath(genpath('cmocean'))

%% generate test scenarios


latlim=[-80 -45];
lonlim=[-65 25];

source_depth=10;

resolution_factor=6
resolution= 2 + ( 10 * (4 ^ resolution_factor) )
[sim_sources.lat,sim_sources.lon] = GridSphere(resolution) ;

% delete grid points ouside study area
ix_inside=sim_sources.lat>latlim(1) & sim_sources.lat<latlim(2) & sim_sources.lon>lonlim(1) & sim_sources.lon<lonlim(2) ;
sim_sources.lat(~ix_inside)=[];
sim_sources.lon(~ix_inside)=[];

% delete grid points on land
[coast]=load('coast');
[Z, R] = vec2mtx(coast.lat, coast.long, ...
    5, [-90 90], [-90 270], 'filled');

val = ltln2val(Z, R, sim_sources.lat,sim_sources.lon);
isOcean = val == 2;

sim_sources.lat(~isOcean)=[];
sim_sources.lon(~isOcean)=[];

load('awi_array_locations.mat')
recorder.lat(8)=[];
recorder.lon(8)=[];
recorder.id=1:numel(recorder.lat);

figure(1)
clf
set(gcf,'color','w')
hold on
worldmap(latlim,lonlim)
[coast]=load('coast');
plotm(coast.lat,coast.long,'-k')
plotm(sim_sources.lat,sim_sources.lon,'.k')

sim_sources.id=1:numel(sim_sources.lat);
sim_sources.id=sim_sources.id';
sim_sources.id(end)
plotm(recorder.lat,recorder.lon,'or')

% get node spacing
for i=1:numel(sim_sources.lat) 
 ix_without_i=1:numel(sim_sources.lat);
ix_without_i(ix_without_i==i)=[];
d_min(i)=min((distance(sim_sources.lat(i), sim_sources.lon(i),sim_sources.lat(ix_without_i), sim_sources.lon(ix_without_i))));
end

node_dist=mean(d_min);
[tl_db]=tl_matrix_geometrical_spreading(sim_sources,recorder,0.150,4000);

%% generate random distribution and save scenario

spacer=.1;
ranlat=[latlim(1):spacer:latlim(2)];
ranlon=[lonlim(1):spacer:lonlim(2)];
[ran.lat,ran.lon]=meshgrid(ranlat,ranlon);

[Z, ran.refvec] = geoloc2grid(ran.lat, ran.lon, ran.lat, spacer);

clear true_sources
true_sources.lat=ran.lat(:)
true_sources.lon=ran.lon(:)
 
% delete grid points on land
[coast]=load('coast');
[Z, R] = vec2mtx(coast.lat, coast.long, ...
    5, [-90 90], [-90 270], 'filled');
val = ltln2val(Z, R, true_sources.lat,true_sources.lon);
isOcean = val == 2;
true_sources.lat(~isOcean)=[];
true_sources.lon(~isOcean)=[];

 [tl.tl_db]=tl_matrix_geometrical_spreading(true_sources,recorder,0.27,4000);

%% scenario loop

n_scenario=1;

for i_scenario=1:n_scenario


ran.x = randnd(-5,size(ran.lat)) ; % Generate 256 x 256 grid with random 1/f (pink) noise
ran.x= ( ran.x-min(ran.x(:)) ) ./ range(ran.x(:)) ;
ran.x(ran.x<.75)=0;
ran.x= ( ran.x-min(ran.x(:)) ) ./ range(ran.x(:)) ;

a=randnd(-1,size(ran.lat));
a= ( a-min(a(:)) ) ./ range(a(:)) ;

ran.x(ran.x>0)=ran.x(ran.x>0) .* a(ran.x>0)
ran.x(ran.x<0)=0;

b=randnd(-2,size(ran.lat));
b= ( b-min(b(:)) ) ./ range(b(:)) ;
ran.x(b<0.6)=0;

true_sources.x=ran.x(:);

% delete grid points on land
[coast]=load('coast');
[Z, R] = vec2mtx(coast.lat, coast.long, ...
    5, [-90 90], [-90 270], 'filled');
val = ltln2val(Z, R, ran.lat,ran.lon);
isOcean = val == 2;

true_sources.x(~isOcean)=[];

sl_mean=180;
call_rate=0.5;
true_sources.sl_p=10.^(sl_mean/20) * call_rate .* true_sources.x ;
true_sources.sl_db=20*log10(true_sources.sl_p);

true_sources.p_true=sum(true_sources.sl_p)

 [p_received,db_received]=received_pressure(true_sources.sl_p,tl.tl_db);
 recorder.p_received=p_received;
 recorder.db_received=db_received;
 
 
ran.p=ran.x*10.^(sl_mean/20) * call_rate ;
ran.db=20*log10(ran.p);
radius = round(node_dist/spacer);
m = fspecial('disk', radius);
ran.smooth = imfilter(ran.p,m,'replicate');
sim_sources.true_p = ltln2val(ran.smooth', ran.refvec, sim_sources.lat,sim_sources.lon);

%%%% save scenario for inversion

save(['scenario_',num2str(i_scenario,'%06.0f')],'ran','sim_sources','recorder','true_sources','tl_db')

%%%%% plot scenario
 
 figure(10)
clf
latlim=[-80 -45];

lonlim=[-65 25];
subplot(211)
hold on
m_proj('lambert','long',lonlim,'lat',latlim);

ran.p=ran.x*10.^(sl_mean/20) * call_rate ;
ran.db=20*log10(ran.p);

m_pcolor(ran.lon,ran.lat,ran.p);
shading flat
m_gshhs_l('patch',[1 1 1]);
m_grid('xlabeldir','end','fontsize',10);
colorbar

plot_bins=100;
plot_spatial= db_received;
plot_range=[min(plot_spatial),max(plot_spatial)];
[~,ind] = histc(plot_spatial,linspace(plot_range(1),plot_range(2),plot_bins));
ind(ind==0)=plot_bins;
cmap=cool(plot_bins);

for i=1:numel(plot_spatial)    
  plot_cmap=cmap(ind(i),:);   
m_plot(recorder.lon(i),recorder.lat(i),'.','markersize',20,'color',plot_cmap)
end

title(['Received db ',num2str(min(db_received)),' - ',num2str(max(db_received))])

subplot(212)
hold on

m_proj('lambert','long',lonlim,'lat',latlim);

m_pcolor(ran.lon,ran.lat,ran.p);
shading flat
m_gshhs_l('patch',[1 1 1]);
m_grid('xlabeldir','end','fontsize',10);

plot_spatial=sim_sources.true_p;
plot_bins=30;
 plot_range=[min(plot_spatial),max(plot_spatial)];
%  plot_range=[0,p_sum_smooth_max];
jet_cmap=cmocean('matter',plot_bins);
[~,ind] = histc(plot_spatial,linspace(plot_range(1),plot_range(2),plot_bins));
 ind(ind==0)=plot_bins;

for j=1:plot_bins
    
plot_cmap=jet_cmap(j,:); 
ix_color=ind==j;

m_plot(sim_sources.lon(ix_color),sim_sources.lat(ix_color),'.','markersize',20,'color',plot_cmap)
end

print(gcf,'-dpng',['scenario_',num2str(i_scenario,'%06.0f')],'-r200')

end

%% estimate sourcre pressure grid using python

% open terminal, navigate to this folder and enter:
% python paraanneal.py

% or

status = system('python paraanneal.py')

%% estimate sourcre pressure grid using matlab

% use the code in sa_parameter_estimation.m
% in and outputs need to be adjusted to the scenario variable names

%% analyse inversion accuracy



sols=dir(['solution_*.mat'])

for i_sol=1:numel(sols)

    
solnume=str2num(sols(i_sol).name(end-9:end-4))
load([sols(i_sol).name]);
load(['scenario',sols(i_sol).name(9:end)]);
ssevalues=sse_iterations(:,end);
if i_sol==1
%%%

% grid of cumulate distance to recs
clear n_recs_in_radius
for isim=1:numel(sim_sources.lat) 
inversion.recdist(i_sol,isim,:)=(deg2km(distance(sim_sources.lat(isim),sim_sources.lon(isim),recorder.lat,recorder.lon)));
inversion.n_recs_in_radius(i_sol,isim)= sum( inversion.recdist(i_sol,isim,:) <= 1000 );    
end
ix=inversion.n_recs_in_radius(i_sol,:)>2;

for isim=1:numel(true_sources.lat) 
inversion.n_recs_in_radius_truesources(isim)=sum(deg2km(distance(true_sources.lat(isim),true_sources.lon(isim),recorder.lat,recorder.lon))  <= 1000  ) ;
inversion.true_sources_cumdist(isim)=sum(deg2km(distance(true_sources.lat(isim),true_sources.lon(isim),recorder.lat,recorder.lon))) / numel(recorder.lat);    
end
ix_ran=inversion.n_recs_in_radius_truesources>=2;
end

ixsourepresent=true_sources.sl_p>0;
inversion.average_dist_to_true_sources(i_sol)=sum(inversion.true_sources_cumdist(ixsourepresent))/sum(ixsourepresent);
inversion.cumulative_dist_to_true_sources(i_sol)=sum(inversion.true_sources_cumdist(ixsourepresent));

%%%%% RL entropy
h=recorder.db_received;
if sum(h)==0
    h=zeros(size(h))+1/numel(h);
else
    h=h/sum(h);
end
en=-sum(h(h>0).*log2(h(h>0)));
b=hist(h,10);
b=b/sum(b);
en_hist=-sum(b(b>0).*log2(b(b>0)));

inversion.rl_entropy(i_sol)=en_hist;
inversion.rl_entropy_spatial(i_sol)=en;

%%%%%%%%% true p  entropy
h=sim_sources.true_p;
if sum(h)==0
    h=zeros(size(h))+1/numel(h);
else
    h=h/sum(h);
end
en=-sum(h(h>0).*log2(h(h>0)));
b=hist(h,10);
b=b/sum(b);
en_hist=-sum(b(b>0).*log2(b(b>0)));
inversion.true_p_entropy(i_sol)=en_hist;
inversion.true_p_entropy_spatial(i_sol)=en;


%%%%%%%%%%% evaluation
[~,ix_best]=sort(likelihood);
ix_best=flip(ix_best);

% estimatied p entropy best
h=est_p(ix_best(1),:);
if sum(h)==0
    h=zeros(size(h))+1/numel(h);
else
    h=h/sum(h);
end
en=-sum(h(h>0).*log2(h(h>0)));
b=hist(h,10);
b=b/sum(b);
en_hist=-sum(b(b>0).*log2(b(b>0)));
inversion.est_p_entropy_best(i_sol)=en_hist;
inversion.est_p_entropy_spatial_best(i_sol)=en;

% estimatied p entropy best3
h=median(est_p(ix_best(1:3),:));
if sum(h)==0
    h=zeros(size(h))+1/numel(h);
else
    h=h/sum(h);
end
en=-sum(h(h>0).*log2(h(h>0)));
b=hist(h,10);
b=b/sum(b);
en_hist=-sum(b(b>0).*log2(b(b>0)));
inversion.est_p_entropy_best3(i_sol)=en_hist;
inversion.est_p_entropy_spatial_best3(i_sol)=en;


%%%%% sse and likelihood
inversion.likelihood_best(i_sol)=likelihood(ix_best(1));
inversion.algorithm_sse_best(i_sol)=sse(ix_best(1));
%best 3
inversion.likelihood_best3(i_sol)=median(likelihood(ix_best(1:3)));
inversion.algorithm_sse_best3(i_sol)=median(ssevalues(ix_best(1:3)));


clear n_recs_in_radius
for isim=1:numel(sim_sources.lat) 
inversion.recdist(i_sol,isim,:)=(deg2km(distance(sim_sources.lat(isim),sim_sources.lon(isim),recorder.lat,recorder.lon)));
inversion.n_recs_in_radius(i_sol,isim)= sum( inversion.recdist(i_sol,isim,:) <= 1000 );    
end
ix=inversion.n_recs_in_radius(i_sol,:)>1;


%%%%%%%%% true p  entropy inside
h=sim_sources.true_p(ix);
if sum(h)==0
    h=zeros(size(h))+1/numel(h);
else
    h=h/sum(h);
end
en=-sum(h(h>0).*log2(h(h>0)));
b=hist(h,10);
b=b/sum(b);
en_hist=-sum(b(b>0).*log2(b(b>0)));
inversion.true_p_entropy_trustzone(i_sol)=en_hist;
inversion.true_p_entropy_spatial_trustzone(i_sol)=en;

% estimatied p entropy best3 inside
h=median(est_p(ix_best(1:3),ix));
if sum(h)==0
    h=zeros(size(h))+1/numel(h);
else
    h=h/sum(h);
end
en=-sum(h(h>0).*log2(h(h>0)));
b=hist(h,10);
b=b/sum(b);
en_hist=-sum(b(b>0).*log2(b(b>0)));
inversion.est_p_entropy_best3_trustzone(i_sol)=en_hist;
inversion.est_p_entropy_spatial_best3_trustzone(i_sol)=en;

%%%%% accuracy best3
x=sim_sources.true_p(ix)';
y=median(est_p(ix_best(1:3),ix));

n_bins=50;
bins=linspace(0,max([x]),n_bins);

[~,ind] = histc(x,bins);
ind(ind==0)=n_bins;
xd=[];
for j=1:numel(bins)
ix_val=ind==j;
xd(ix_val)=bins(j);
end

 bins=linspace(0,max([x]),n_bins); % testing -- bad idea

[~,ind] = histc(y,bins);
ind(ind==0)=n_bins;
yd=[];
for j=1:numel(bins)
ix_val=ind==j;
yd(ix_val)=bins(j);
end

intersect=sum(xd==yd);
inversion.accuracy_norm(i_sol)=intersect / numel(x);

%binary accuracy best3
xb=x;
yb=y;
xb(xb>0)=1;
yb(yb>0)=1;

intersect=sum(xb==yb);
inversion.accuracy_binary(i_sol)= intersect / numel(x);

%entiregrid
x=sim_sources.true_p';
y=median(est_p(ix_best(1:3),:));

xb=x;
yb=y;
xb(xb>0)=1;
yb(yb>0)=1;

intersect=sum(xb==yb);
inversion.accuracy_binary_grid(i_sol)= intersect / numel(x);


n_bins=50;
bins=linspace(0,max([x]),n_bins);

[~,ind] = histc(x,bins);
ind(ind==0)=n_bins;
xd=[];
for j=1:numel(bins)
ix_val=ind==j;
xd(ix_val)=bins(j);
end

[~,ind] = histc(y,bins);
ind(ind==0)=n_bins;
yd=[];
for j=1:numel(bins)
ix_val=ind==j;
yd(ix_val)=bins(j);
end

intersect=sum(xd==yd);
inversion.accuracy_norm_grid(i_sol)=intersect / numel(x);

%%%%%
%accuracy of best sol

x=sim_sources.true_p(ix)';
y=est_p(ix_best(1),ix);

n_bins=50;
bins=linspace(0,max([x]),n_bins);

[~,ind] = histc(x,bins);
ind(ind==0)=n_bins;
xd=[];
for j=1:numel(bins)
ix_val=ind==j;
xd(ix_val)=bins(j);
end

[~,ind] = histc(y,bins);
ind(ind==0)=n_bins;
yd=[];
for j=1:numel(bins)
ix_val=ind==j;
yd(ix_val)=bins(j);
end

intersect=sum(xd==yd);
inversion.accuracy_norm_best(i_sol)=intersect / numel(x);

%binary accuracy best
xb=x;
yb=y;
xb(xb>0)=1;
yb(yb>0)=1;

intersect=sum(xb==yb);
inversion.accuracy_binary_best(i_sol)= intersect / numel(x);

% % %%%%%%%%%%%%%%%%%%%%%%%


inversion.rl(i_sol,:)=recorder.db_received;
inversion.rec_lat(i_sol,:)=recorder.lat;
inversion.rec_lon(i_sol,:)=recorder.lon;

inversion.t(i_sol,:)=sim_sources.true_p;
inversion.sim_lat(i_sol,:)=sim_sources.lat;
inversion.sim_lon(i_sol,:)=sim_sources.lon;

inversion.s_best3(i_sol,:)=median(est_p(ix_best(1:3),:))';
inversion.s_best(i_sol,:)=est_p(ix_best(1),:)';

inversion.ran_p(i_sol,:,:)=ran.p;
inversion.ran_lat(i_sol,:,:)=ran.lat;
inversion.ran_lon(i_sol,:,:)=ran.lon;

%%%%%%%%%%%%%%%%%
inversion.true_p_sum(i_sol)=true_sources.p_true;
inversion.est_p_sum(i_sol)=sum(median(est_p(ix_best(1:3),:)));

inversion.true_p_sum_trustzone(i_sol)=sum(true_sources.sl_p(ix_ran));
inversion.true_node_p_sum_trustzone(i_sol)=sum(sim_sources.true_p(ix));
inversion.est_p_sum_trustzone(i_sol)=sum(median(est_p(ix_best(1:3),ix)));

end

