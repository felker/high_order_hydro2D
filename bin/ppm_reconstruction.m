%ppm_reconstruction.m
%plot the exact piecewise parabolas of the reconstructed state
close all; 
title_str = sprintf('PPM reconstructed states, cycle=0'); %have to edit manually for now
var = dlmread('ppm.txt',',',0,0);
xrange_view = 22:26; %4:259;  %full view excludes ghost cell ie+1
xrange_full = 4:260; %currently, cout in van_leer2.cpp goes from is-1:ie+1, 
%3:260 for NGHOST=4 nx1=256, is=4, ie=259
%STILL DONT KNOW WHY i output RIEMANN INTERFACE STATES FOR IS-1.
%they arent even set inside PLM.cpp, ppm.cpp
nsamples = 20; %invisible points inside cell boundaries used to create continuous parabola
ylab_precision = 12;

prim_char = char('density \rho','v_x','v_y','pressure P');
prim_str = cellstr(prim_char); %this becomes a cell array of strings (since they have different lengths!)
color_char = char('m','r','k','g');

VAR_INDEX = 3; 
%new order of variables in this CSV output is w, wl, wr ...
%where wl,wr are riemann states at i-1/2 boundary, NOT parabolic
%interpolation limits
%16 columns 
if (VAR_INDEX==4) %skip vz to get pressure
    w = var(xrange_full-2,14);
    wl = var(xrange_full-2,15); 
    wr = var(xrange_full-2,16);        
else
    w = var(xrange_full-2,(VAR_INDEX-1)*3 +1);
    wl = var(xrange_full-2,(VAR_INDEX-1)*3+2); 
    wr = var(xrange_full-2,(VAR_INDEX-1)*3+3); 
end
%now, w variables go from 1 to 257, for the one ghost cell after the x1
%boundary
%convert to al, ar for all real cells, is:ie, 4:259
al = wr(1:256); 
ar = wl(2:257); 
a = w(1:256); 
delta_a = ar - al; 
a_6 = 6*(a - 0.5*(al+ar)); 

figure; 
hold on;
%xrange_view(i) is the ATHENA++ cell number, but the corresponding i in the
%loop below refers to al(i) at i-1/2, ar(i) at i+1/2, i.e. I choose to cell
%center the index (could have also centterd index at interface
for i=1:size(xrange_view')
    %here, "xi" is not physical x position, but athena++ i index
    %dxi = 1
    xi = linspace(xrange_view(i)-0.5,xrange_view(i)+0.5,nsamples);  
    %use this variable to index the a variables
    index = xrange_view(i) - 3; 
    x=(xi-xi(1));
    a_reconstruct = al(index) + x.*(delta_a(index) + a_6(index)*(1-x)); 
    plot(xi,a_reconstruct,'-k'); 
    %plot interface points
    plot(xi(1),al(index),'*k'); 
    plot(xi(1)+1,ar(index),'*k'); 
    %plot cell average
    plot(xi(1)+0.5,a(index),'*b'); 
    %plot(xrange_view(i)-0.5,rho(xrange_view(i)-2),'*r'); %interface values rho go from 3:261
end
xlabel('ATHENA++ i-index'); 
legend('Parabolic interpolant','a_{-/L}','a_{+/R}','<a>_j');
yt=get(gca,'YTick');
ylab=num2str(yt(:), ylab_precision);
set(gca,'YTickLabel',ylab);