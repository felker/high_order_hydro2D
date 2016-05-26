%reconstructed_states_reader.m
%script for reading and plotting the CSV file that is produced from my 
%van_leer2.cpp debugging cout statements when stripped down to a single
%timestep. See formatted_table_reader.m for some comments about these
%variables

%ohhhhh the string arrays only get the first letter of each word if you
%access it like a matrix-- use cellstr() to 
%Convert to cell array of character vectors

% change "figure(3)" to something more descriptive
% make figures tile instead of on top of each other

% 
% Specifying Figure Size and Screen Location
% To create a figure window that is one quarter the size of your screen and is positioned in the upper left corner, use the root object's ScreenSize property to determine the size. The ScreenSize property value is a four-element vector: [left bottom width height].
% 
% scrsz = get(groot,'ScreenSize');
% figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])
% To position the full figure window including the menu bar, title bar, tool bars, and outer edges, use the OuterPosition property in the same manner.
% 
% Specifying the Figure Window Title
% You can add your own title to a figure by setting the Name property and turning off the NumberTitle property:
% 
% figure('Name','Simulation Plot Window','NumberTitle','off')


close all;
hold off; 

%make variable switching automatic

%var = plm; % name of variable I imported
title_str = sprintf('PPM reconstructed states, cycle=0'); %have to edit manually for now
var = dlmread('ppm.txt',',',0,0);
xrange = 4:259; %22:26; %currently, cout goes from is-1:ie+1, 3:260 for NGHOST=4
ylab_precision = 8;

prim_char = char('density \rho','v_x','v_y','pressure P');
prim_str = cellstr(prim_char); %this becomes a cell array of strings (since they have different lengths!)
color_char = char('m','r','k','g');

%Use this section to print out all 4 primitive variables
for i=1:4
    ax(i) = subplot(2,2,i); 
    if (i==4) %skip vz to get pressure
        plot(ax(i),xrange,var(xrange-2,9),'--o','Color',color_char(i));
        hold on;
        plot(ax(i),xrange,var(xrange-2,10),':o','Color',color_char(i));        
    else
        plot(ax(i),xrange,var(xrange-2,(i-1)*2 +1),'--o','Color',color_char(i));
        hold on;
        plot(ax(i),xrange,var(xrange-2,(i-1)*2+2),':o','Color',color_char(i));
    end
    title(ax(i),prim_str(i),'Interpreter','tex');
    legend(ax(i), 'L','R'); 
    yt=get(ax(i),'YTick');
    ylab=num2str(yt(:), ylab_precision);
    set(ax(i),'YTickLabel',ylab);
    xlabel('ATHENA++ i-index');     
end

%plot global title
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,title_str,'HorizontalAlignment' ,'center','VerticalAlignment', 'top')

%link x axes so I can zoom simultaneously
linkaxes(ax(1:4), 'x');
%END SECTION

%use this section to only examine one variable in detail
VAR_INDEX = 2; 
figure; 
if (VAR_INDEX==4) %skip vz to get pressure
    plot(xrange,var(xrange-2,9),'--o','Color',color_char(VAR_INDEX));
    hold on;
    plot(xrange,var(xrange-2,10),':o','Color',color_char(VAR_INDEX));        
else
    plot(xrange,var(xrange-2,(VAR_INDEX-1)*2 +1),'--o','Color',color_char(VAR_INDEX));
    hold on;
    plot(xrange,var(xrange-2,(VAR_INDEX-1)*2+2),':o','Color',color_char(VAR_INDEX));
end
title(prim_str(VAR_INDEX),'Interpreter','tex');
legend('L','R'); 
yt=get(gca,'YTick');
ylab=num2str(yt(:), ylab_precision);
set(gca,'YTickLabel',ylab);
xlabel('ATHENA++ i-index'); 


%%%%%IMPORTANT%%%%%%
%quit here after first run through
ppm = var; 
return;

%USE THIS SECTION TO COMPARE PLM AND PPM RECONSTRUCTIONS FOR ONE VARIABLE
%run the above code twice, rename var appropriately
var1 = ppm; 
var2 = plm; 
figure; 
if (VAR_INDEX==4) %skip vz to get pressure
    %temporarily change the color for PPM
    plot(xrange,var1(xrange-2,9),'--o','Color','b');
    hold on;
    plot(xrange,var1(xrange-2,10),':o','Color','b');  
    plot(xrange,var2(xrange-2,9),'-o','Color',color_char(VAR_INDEX));
    plot(xrange,var2(xrange-2,10),'-.o','Color',color_char(VAR_INDEX));
else
%temporarily change the color for PPM 
%     plot(xrange,circshift(var1(xrange-2,(VAR_INDEX-1)*2+1),0),'--o','Color','b');
%     hold on;
%     plot(xrange,circshift(var1(xrange-2,(VAR_INDEX-1)*2+2),0),':o','Color','b');
    plot(xrange,var1(xrange-2,(VAR_INDEX-1)*2+1),'--o','Color','b');
    hold on;
    plot(xrange,var1(xrange-2,(VAR_INDEX-1)*2+2),':o','Color','b');
    plot(xrange,var2(xrange-2,(VAR_INDEX-1)*2+1),'-o','Color',color_char(VAR_INDEX));
    plot(xrange,var2(xrange-2,(VAR_INDEX-1)*2+2),'-.o','Color',color_char(VAR_INDEX));
end
title(prim_str(VAR_INDEX),'Interpreter','tex');
legend('PPM L','PPM R','PLM L','PLM R'); 
yt=get(gca,'YTick');
ylab=num2str(yt(:), ylab_precision);
set(gca,'YTickLabel',ylab);
xlabel('ATHENA++ i-index'); 

% %Plot separation of L/R states. Very hacked-up. Need to run twice, changing
% %the file above for PLM and PPM, and 
figure;
if (VAR_INDEX==4) %skip vz to get pressure
    ppm_state_distance = abs(var1(:,9) - var1(:,10));
    plm_state_distance = abs(var2(:,9) - var2(:,10));
else
    ppm_state_distance = abs(var1(:,(VAR_INDEX-1)*2 +1) - var1(:,(VAR_INDEX-1)*2+2));
    plm_state_distance = abs(var2(:,(VAR_INDEX-1)*2 +1) - var2(:,(VAR_INDEX-1)*2+2));
end
% future fix: make this above difference span all prim variables
semilogy(xrange,ppm_state_distance(xrange-2),'-o','Color',color_char(VAR_INDEX))
hold on; 
semilogy(xrange,plm_state_distance(xrange-2),':o','Color',color_char(VAR_INDEX));
legend('PPM','PLM'); 
xlabel('ATHENA++ i-index'); 
title(sprintf('| %s_L - %s_R |',prim_str{VAR_INDEX},prim_str{VAR_INDEX})); 


% %Plot separation of PPM_L and PLM_L states, in addition to R
figure;
if (VAR_INDEX==4) %skip vz to get pressure
    R_state_distance = abs(var2(:,10) - var1(:,10));
    L_state_distance = abs(var2(:,9) - var1(:,9));
else
    R_state_distance = abs(var2(:,(VAR_INDEX-1)*2 +2) - var1(:,(VAR_INDEX-1)*2+2));
    L_state_distance = abs(var2(:,(VAR_INDEX-1)*2 +1) - var1(:,(VAR_INDEX-1)*2+1));
end
% future fix: make this above difference span all prim variables
semilogy(xrange,L_state_distance(xrange-2),'-o','Color',color_char(VAR_INDEX))
hold on; 
semilogy(xrange,R_state_distance(xrange-2),':o','Color',color_char(VAR_INDEX));
legend('L state','R state'); 
xlabel('ATHENA++ i-index'); 
title(sprintf('| PPM %s_L - PLM %s_L |',prim_str{VAR_INDEX},prim_str{VAR_INDEX})); 