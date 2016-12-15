function [d]=Homer2SPMdata(d)

fn = fullfile(d.root,d.babydir,d.raw)

% import the raw data from Homer
% the hd matrix contains all the info about raw data
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

%% channel definition
%% source, detector, channel in the order from hd.SD.MeasList but with channel numbers added
d.chlab = [3   7   15
    3   8   16
    4   8   19
    4   10  20
    4   11  22
    5   1   2
    5   3   3
    6   1   1
    6   3   4
    6   4   5
    8   7   14
    8   8   17
    8   10  18
    9   11  25 
    9   15  26
    10  3   6
    10  4   7
    10  6   9
    13  4   8
    13  6   10
    13  16  11
    15  10  21
    15  11  23
    15  15  24
    16  6   12
    16  16  13];

%% d.chlab from hd
% d.chlab(:,1) = hd.SD.MeasList(1:26,1)
% d.chlab(:,2) = hd.SD.MeasList(1:26,2)
% 3rd col??
cbad = (zeros(d.Nch,1));

%importing the channels marked in yellow (Homer)
if ~isempty(hd.procResult.SD)
    tmp = hd.procResult.SD.MeasListAct; 
    tmp = reshape(tmp,d.Nch,2);
    
    cbad1 = all(tmp>0.5,2); 
end

%importing the channels marked in pink (KP script) %%hd.SD.MeasListAct
%(hd.procInput.SD.MeasListAct)
if ~isempty(hd.SD)
    tmp2 = hd.SD.MeasListAct; 
    tmp2 = reshape(tmp2,d.Nch,2);
    
    cbad2 = all(tmp2>0.5,2); 
end

%creating the cbad variable with bad channels (both yellow and pink ones)
cbad = cbad1 & cbad2;

%keyboard 
%???
[i,reorder] = sort(d.chlab(:,3));

%importing oxy and deoxy data from Homer (not sure about reorder) ??
d.nirs_data.oxyData = squeeze(hd.procResult.dc(:,1,reorder));  
d.nirs_data.dxyData = squeeze(hd.procResult.dc(:,2,reorder)); 

%% store all the info in the d matrix
d.nirs_data.oxyData=d.nirs_data.oxyData *1000000;
d.nirs_data.dxyData=d.nirs_data.dxyData *1000000;

%% reorder cbad to match channels
%cbad = cbad(reorder);

d.nirs_data.nchn = d.Nch;
d.nirs_data.fs = 10;  %% frequency
d.nirs_data.wavelength = [770 850];
d.nirs_data.distance = ones(1,d.nirs_data.nchn)*2.0; % change!!! 
d.nirs_data.DPF = 5.1300;

d.nirs_data.DPF_correction= 'Charite correction';
d.rate = 10;% KB
d.cbad=logical(cbad);
d.cname = [1:d.Nch];


%keyboard

figure(1), clf % figure of the channels structure
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

%%save the picture of bad channels in the subject folder
pname = fullfile(d.root,d.babydir,'bad_ch.tif');
print('-dtiff','-zbuffer',pname)

cbad = cbad(reorder);
d.cbad=logical(cbad);
