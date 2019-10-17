% Voice Analysis
%
% Jorge C. Lucero
% 
% Voice analysis based on: H. Herzel, D. Berry, I. R. Titze, & M. Saleh, 
% “Analysis of vocal disorders with methods from nonlinear dynamics”, 
% Journal of Speech and Hearing Research 37, 1008-1019, 1994.

clear all

FFI = 220;  %This must be slightly higher than the voice F0

% Input voice record (wav)

[file,path] = uigetfile('*.wav');
if isequal(file,0)
   disp('Cancel');
else
   disp(['File selected: ', fullfile(path,file)]);
end

%[y, Fs] = wavread(file);
[y, Fs] = audioread(file); 
len = length(y);

% Low pass filter

Ws = 2*FFI/Fs;
[b, a] = butter(5, Ws);
yyfilt = filtfilt(b, a, y);

% Analyze F0

ind = 0;
for n = 1:len-1
    if(yyfilt(n) <= 0 & yyfilt(n + 1) >= 0)
        ind = ind + 1;
        Pp(ind) = n;
        Pr(ind) = n - yyfilt(n)/(yyfilt(n + 1) - yyfilt(n));
    end
end

Pp(1) = [];
Pr(1) = [];
ind = ind-1;

for n = 1:ind - 1
    f0(n) = Fs/(Pr(n + 1) - Pr(n));
end

meanf0 = mean(f0);
stdf0 = std(f0);
maxf0 = max(f0);
minf0 = min(f0);

if maxf0 >= FFI
    warning(['Low-pass cutoff frequency ', num2str(FFI),...
          '(Hz) is not larger than max F0 ',num2str(maxf0),...
          '(Hz)'])
end
          
pert0 = (f0(1:ind - 1) - meanf0)/meanf0;
pert1 = pert0(2:ind - 1) - pert0(1:ind - 2);
pert2 = pert1(2:ind - 2) - pert1(1:ind - 3);

mr1 = mean(abs(pert1));
mr2 = mean(abs(pert2));

halfind = floor((ind - 1)/2);
for i = 0:halfind;
    f0corr(i + 1) = sum(pert0(1:halfind).*pert0(1 + i:halfind + i));
end
f0corr = f0corr./f0corr(1);


% Analyze amplitude 

for n = 1:ind - 1
    amp(n) = max(yyfilt(Pp(n):Pp(n + 1))) - min(yyfilt(Pp(n):Pp(n + 1)));
end

meanamp = mean(amp);
stdamp = std(amp);
maxamp = max(amp);
minamp = min(amp);

apert0 = (amp(1:ind - 1) - meanamp)/meanamp;
apert1 = apert0(2:ind - 1) - apert0(1:ind - 2);
apert2 = apert1(2:ind - 2) - apert1(1:ind - 3);

amr1 = mean(abs(apert1));
amr2 = mean(abs(apert2));

for i = 0:halfind;
    ampcorr(i + 1) = sum(apert0(1:halfind).*apert0(1 + i:halfind + i));
end
ampcorr = ampcorr./ampcorr(1);
 
% Plot results 

h1 = figure(1);
clf;

hs1 = subplot(2, 1, 1);
plot((0:len-1)/Fs, y);
set(hs1, 'Fontsize', 10);
title(['Acoustics - File: ',file]);
xlabel('Time (s)');

hs2 = subplot(2, 1, 2);
spectrogram(y, round(.1*Fs), round(.08*Fs), round(.1*Fs), Fs, 'yaxis')
load cmap
colormap(cmap);
colorbar('off');
grid off
xlabel('Time (s)');

linkaxes([hs1 hs2], 'x');

h2=figure(2);
clf;
pos = get(h2,'Position');
set(h2,'Position',[(pos(1) + 10) (pos(2) - 110) 560 520]);

hs1=subplot(3, 1, 1);
plot(f0)
v = axis;
v(2) = ind - 1;
axis(v);
set(hs1, 'Fontsize', 10);
title('F0 contour');
xlabel('No. cycles');
ylabel('Hz');
grid on

hs2 = subplot(3, 1, 2);
plot(pert2*100)
v = axis;
v(2) = ind - 1;
axis(v);
set(hs2, 'Fontsize', 10);
title('2nd perturbation');
xlabel('No. cycles');
ylabel('%');
grid on;

hs3 = subplot(3, 2, 5);
plot(f0corr);
v = axis;
v(2) = length(f0corr);
v(3) = -1;
v(4) = 1;
axis(v);
set(hs3, 'Fontsize', 10);
title('Autocorrelation');
xlabel('No. cycles');
grid on

hs4 = subplot(3, 2, 6);
set(hs4, 'Visible', 'off');
set(hs4, 'Fontsize', 10);
text(0, 1, ['File = ', file], 'Fontsize', 10);
text(0, .8, ['Mean F0 = ', num2str(meanf0), ' Hz'], 'Fontsize', 10);
text(0, .65, ['Max  F0 = ', num2str(maxf0), ' Hz'], 'Fontsize', 10);
text(0, .5, ['Min  F0 = ', num2str(minf0), ' Hz'], 'Fontsize', 10);
text(0, .35, ['Std  F0 = ', num2str(stdf0), ' Hz'], 'Fontsize', 10);
text(0, .2, ['MR1     = ', num2str(mr1*100), '%'], 'Fontsize', 10);
text(0, .05, ['MR2     = ', num2str(mr2*100), '%'], 'Fontsize', 10);

h3=figure(3);
clf;
pos = get(h3,'Position');
set(h3,'Position', [(pos(1) + 20) (pos(2) - 120) 560 520]);

hs1 = subplot(3, 1, 1);
plot(amp)
v = axis;
v(2) = ind - 1;
axis(v);
set(hs1, 'Fontsize', 10);
title('Amplitude contour');
xlabel('No. cycles');
grid on

hs2 = subplot(3, 1, 2);
plot(apert2*100)
v = axis;
v(2) = ind - 1;
axis(v);
set(hs2, 'Fontsize', 10);
title('2nd perturbation');
xlabel('No. cycles');
ylabel('%');
grid on;

hs3 = subplot(3, 2, 5);
plot(ampcorr);
v = axis;
v(2) = length(ampcorr);
v(3) = -1;
v(4) = 1;
axis(v);
set(hs3, 'Fontsize', 10);
title('Autocorrelation');
xlabel('No. cycles');
grid on

hs4 = subplot(3, 2, 6);
set(hs4, 'Visible', 'off');
set(hs4, 'Fontsize', 10);
text(0, 1, ['File = ', file], 'Fontsize', 10);
text(0, .8, ['Mean Amp = ', num2str(meanamp)], 'Fontsize', 10);
text(0, .65, ['Max  Amp = ', num2str(maxamp)], 'Fontsize', 10);
text(0, .5, ['Min  Amp = ', num2str(minamp)], 'Fontsize', 10);
text(0, .35, ['Std  Amp = ', num2str(stdamp)], 'Fontsize', 10);
text(0, .2, ['MR1     = ', num2str(amr1*100), '%'],'Fontsize', 10);
text(0, .05, ['MR2     = ', num2str(amr2*100), '%'],'Fontsize', 10);
