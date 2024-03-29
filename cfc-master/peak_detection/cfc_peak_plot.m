function cfc_peak_plot(obj, outpath, max_freq)
%%
% Create a plot from the output struct from sgolay_peaks
%
% obj: peaks struct
%   the output from a call to cfc_peak_detect
% outpath: string
%   a filepath to save out a copy of the figure
% maxfreq: double
%   the maximum frequency you want to display

if nargin<3
    max_freq = obj.freq_vect(end-1);
end

if nargin<2
    outpath =false;
end

max_freq_idx = length(obj.freq_vect(obj.freq_vect < max_freq));

figure('Position',[50 50 500 800]);
subplot(311)
plot(obj.freq_vect(1:max_freq_idx),obj.psd(1:max_freq_idx))
hold on
plot(obj.freq_vect(1:max_freq_idx),obj.smo_psd(1:max_freq_idx),'r')

for ipk = 1:length(obj.peak)
    plot(obj.peak{ipk}.peak_freq,obj.peak{ipk}.peak_amp,'og');
end
title('Frequency Spectrum')
legend({'Original','Smoothed','Peaks'})
ylabel('Amplitude')
grid on;

subplot(312)
plot(obj.freq_vect(2:max_freq_idx),diff(obj.smo_psd(1:max_freq_idx)))
hold on
grid on;

for ipk = 1:length(obj.peak)
    if obj.peak{ipk}.p < (0.05 / length(obj.peak));
        plot(obj.peak{ipk}.X(1,:),obj.peak{ipk}.Y_hat,'g','linewidth',2);
    end
end
title('Differential of spectrum')
xlabel('Frequency (Hz)')
legend({'Original','Linear fits'})

%% Fancy table
table = uicontrol('style','text','Position',[25 50 450,200],'FontName','Courier New');
header = sprintf('%-8s%-10s%-12s%-12s%-8s%-8s%-8s%-8s','Peak','Freq(Hz)','Amp','Beta','R2','SE','T','p');
table_lines = repmat('-',1,length(header));
table_content = strcat(table_lines,'\n',header,'\n',table_lines);
for idx = 1:length(obj.peak)
    table_content = strcat(table_content,'\n',...
        sprintf('%-8d%-10.2f%-12.2f%-12.2f%-8.2f%-8.2f%-8.2f%-8.2f\n',...
        idx,...
        obj.peak{idx}.peak_freq,...
        obj.peak{idx}.peak_amp,...
        obj.peak{idx}.coeff(1),...
        obj.peak{idx}.r2,...
        obj.peak{idx}.SE,...
        obj.peak{idx}.t,...
        obj.peak{idx}.p));
end
table_content = strcat(table_content,'\n',table_lines);
set(table,'String',sprintf(table_content));


if outpath ~= false
    savefig(outpath);
    print(gcf, outpath, '-deps');
end

