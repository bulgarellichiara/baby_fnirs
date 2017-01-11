function [d]=Homer2SPMdata(d)
%% last version of fload_homer (11-01-17 -- CB)
% right reording of bad channels, output figure with bad channels marked in,
% all the info are read in from Homer, no need to specify anything!
% !! please create an excel file with the list of the channel named 'ch.xlsx' and
% save it in the d.root folder

%%%%%% 

fn = fullfile(d.root,d.babydir,d.raw)

% import the all the raw data from Homer in the hd matrix
hd = importdata(fn)

%% get bad trials
badt = cat(1,hd.userdata.data{:,1});

for i=1:length(badt)
    mark = hd.userdata.data{i,2};
    if(length(mark)>0)
        badt(i,2) = 1;
    else
        badt(i,2) = 0;
    end
end
    
d.badt = badt;  %% store bad trials in the d matrix

et = hd.s;

ett = zeros(length(et),1);
%% add up events
for i=1:size(et,2)
    ett = ett+et(:,i).*i;
end

d.ev = ett;  %% very basic event coding by Homer column
d.condnames = hd.CondNames;  %% read from Homer
  
d.Scrpos = hd.SD.SrcPos; %info about sources
d.Detpos = hd.SD.DetPos; %info about detectors
d.Nch = length(hd.d(1,:))/2; % number of channels detected automatically, no need to be changed!


%% get sources, detectors, channels in the d.chlab from hd
d.chlab(:,1) = hd.SD.MeasList(1:d.Nch,1);
d.chlab(:,2) = hd.SD.MeasList(1:d.Nch,2);
ch = dir(fullfile(d.root,'ch.xlsx'));
d.chlab(:,3)=xlsread(ch.name);
d.chlab
%check the d.chlab printed out in matlab 

cbad = (zeros(d.Nch,1));

%importing the channels marked in yellow (from Homer)
if ~isempty(hd.procResult.SD)
    tmp = hd.procResult.SD.MeasListAct; 
    tmp = reshape(tmp,d.Nch,2);
    
    cbad1 = all(tmp>0.5,2); 
end

%importing the channels marked in pink (from KP script)
if ~isempty(hd.SD)
    tmp2 = hd.SD.MeasListAct; 
    tmp2 = reshape(tmp2,d.Nch,2);
    
    cbad2 = all(tmp2>0.5,2); 
end

%creating the cbad variable with bad channels (both yellow and pink ones)
cbad = cbad1 & cbad2;
[i,reorder] = sort(d.chlab(:,3));

%importing oxy and deoxy data from Homer 
d.nirs_data.oxyData = squeeze(hd.procResult.dc(:,1,reorder));  
d.nirs_data.dxyData = squeeze(hd.procResult.dc(:,2,reorder)); 

% store all the info about the data in the d matrix
d.nirs_data.oxyData=d.nirs_data.oxyData *1000000;
d.nirs_data.dxyData=d.nirs_data.dxyData *1000000;

d.nirs_data.DPF_correction= 'Charite correction';
d.rate = 10;% KB
d.cbad=logical(cbad);
d.cname = [1:d.Nch];

%% figure of the channels structure with bad channels marked in black
figure(1), clf 
for i=1:length(d.chlab)
    ss = d.Scrpos(d.chlab(i,1),:);
    dd = d.Detpos(d.chlab(i,2),:);
    plot([ss(1),dd(1)], [ss(2),dd(2)], 'b-','Color',[0.2,0.2,0.8])
    hold on
    h=plot( mean([ss(1),dd(1)]) , mean([ss(2),dd(2)]) , 'bo','MarkerSize',15,'MarkerFaceColor',[1,1,1],'Color',[0.2,0.2,0.8]);
    text( mean([ss(1),dd(1)]) , mean([ss(2),dd(2)]) , num2str(d.chlab(i,3)),'HorizontalAlignment','Center');

    if(d.cbad(i)==0)
        h=plot( mean([ss(1),dd(1)]) , mean([ss(2),dd(2)]) , 'bo','MarkerSize',15,'MarkerFaceColor',[0,0,0],'Color',[0,0,0]);
        text( mean([ss(1),dd(1)]) , mean([ss(2),dd(2)]) , num2str(d.chlab(i,3)),'HorizontalAlignment','Center','Color',[1,1,1]);
    end        
end
axis off

% sources are marked in red, while detectors are marked in green
plot(d.Scrpos(:,1),d.Scrpos(:,2),'ro','MarkerFaceColor',[1,0,0],'MarkerSize',8);
hold on
for i=1:length(d.Scrpos)
    text(d.Scrpos(i,1)+1,d.Scrpos(i,2)+1,num2str(i))
end
plot(d.Detpos(:,1),d.Detpos(:,2),'go','MarkerFaceColor',[0,1,0],'MarkerSize',8);
hold on
for i=1:length(d.Detpos)
    text(d.Detpos(i,1)+1,d.Detpos(i,2)+1,num2str(i))
end

%%save the picture of channels structures and bad channels in the subject folder
pname = fullfile(d.root,d.babydir,'bad_ch.tif');
print('-dtiff','-zbuffer',pname)

%saving bad channels in the right order in the d matrix
cbad = cbad(reorder);
d.cbad=logical(cbad);