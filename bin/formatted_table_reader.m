%formatted_table_reader.m
%tool for reading in .tab file from ATHENA++, sliced in 1D

%observations/ todo:
% cant zoom in on subplots!!
% subscripted subplot titles wont appear
% fix window size to be larger?
% upgrade to load many timeslices?
% make general for other slices
% fix the cycle=0 header bug?

close all 
%clear all
format LONGG

%USER OPTIONS
cycle=1; %integer cycle  >= 0 at which you want to probe the lineout
digits=5; %length field of format specifier of .tab file numbers. 
ylab_precision = 8; %integer number of digits after decimal to show on yaxis labels
%e.g. LinWave.block0.out2.00000.tab. Not sure how this is set in ATHENA++? worth learning
xrange = 4:259; %10:35;  %22:25;%%specify ATHENA++ i-indices
%END

%recall, real cell index starts at 4 currently in ATHENA++
NGHOST=4; 
xrange = xrange - (NGHOST-1); %since MATLAB is 1-indexed, C++ is 0 indexed

title_str = sprintf('PPM cycle = %d',cycle);
format_str = sprintf('%%0%d.0f',digits); 
%cool! you can nest sprintf() calls to make the format specifier a variable
file_str = sprintf(sprintf('LinWave.block0.out2.%s.tab',format_str),cycle) 
%two row offset otherwise
%e.g. 
%# Athena++ data at time=0  cycle=0  variables=prim
%# Slice at x3=0  (k-ks)=0
%# Slice at x2=0.187795  (j-js)=21
%v.s.
%# Athena++ data at time=0.0331916  cycle=19  variables=prim
%# Slice at x2=0.187795  (j-js)=21
if (cycle ==0) %can probably simplify by eliminating x3 slice option in output
    row_offset=3; 
else
    row_offset=2;
end
LinWave = dlmread(file_str,' ',row_offset,0); %3 row offset for initial condition cycle=0

%automatically plot all 4 primitive 2D variables
figure
prim_char = char('','','density \rho','pressure P','vx','vy','vz');
prim_str = cellstr(prim_char);
color_char = char('','','m','g','k','r','');
for i=3:6
    ax(i-2) = subplot(2,2,i-2); 
    plot(ax(i-2),LinWave(xrange,1),LinWave(xrange,i),'-o','Color',color_char(i));
    title(ax(i-2),prim_str(i),'Interpreter','tex');
    yt=get(ax(i-2),'YTick');
    ylab=num2str(yt(:), ylab_precision);
    set(ax(i-2),'YTickLabel',ylab);
    xlabel('ATHENA++ i-index'); 
end

%plot global title
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,title_str,'HorizontalAlignment' ,'center','VerticalAlignment', 'top')

%link x axes so I can zoom simultaneously
linkaxes(ax(1:4), 'x');

%%%%%IMPORTANT%%%%%%
%quit here after first run through (until both ppm and plm are assigned)
ppm = LinWave; 

%USE THIS SECTION TO COMPARE PLM AND PPM SOLUTIONS FOR ONE VARIABLE
%run the above code twice, rename plm, ppm as LinWave appropriately

VAR_INDEX=5; %between 3 and 6 here, not 1:4 as in other reader. diff order
var1 = ppm; 
var2 = plm; 
figure; 

plot(var1(xrange,1),var1(xrange,VAR_INDEX),'--o','Color',color_char(VAR_INDEX));
hold on;
plot(var2(xrange,1),var2(xrange,VAR_INDEX),':o','Color',color_char(VAR_INDEX));

title(prim_str(VAR_INDEX),'Interpreter','tex');
legend('PPM','PLM'); 
yt=get(gca,'YTick');
ylab=num2str(yt(:), ylab_precision);
set(gca,'YTickLabel',ylab);
xlabel('ATHENA++ i-index'); 

%plot difference between plm and ppm solutions
figure; 
semilogy(var2(xrange,1),abs(var2(xrange,VAR_INDEX)-var1(xrange,VAR_INDEX)) ,'-o','Color',color_char(VAR_INDEX));
xlabel('ATHENA++ i-index'); 
title(sprintf('| %s_{PLM} - %s_{PPM} | cycle= %d solution',prim_str{VAR_INDEX},prim_str{VAR_INDEX},cycle)); 
