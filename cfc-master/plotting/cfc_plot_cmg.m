function cfc_plot_cmg(cmg,metric_name,spectrum,freq_vect,log_high,col_levels)

if nargin < 2 || isempty(metric_name)
    if isfield(cmg,'mi_norm')
        metric_name = 'mi_norm';
    elseif isfield(cmg,'mi')
        metric_name = 'mi';
    elseif isfield(cmg,'esc')
        metric_name = 'esc';
    elseif isfield(cmg,'nesc')
        metric_name = 'nesc';   
   elseif isfield(cmg,'glm')
        metric_name = 'glm';  
      elseif isfield(cmg,'aec')
        metric_name = 'aec';      
           elseif isfield(cmg,'plv')
        metric_name = 'plv';      
        
    else
        error('Please specify which metric to plot!')
    end
end

metric = cmg.(metric_name);

if ndims(metric) >= 3
    metric = squeeze(mean(metric,4));
end

if nargin <= 2

    %% Just plot the cmg without any fancy spectrum sidebars
    col_levels = 0:max(max(metric))/24:max(max(metric));
    figure;
    contourf(squeeze(mean(cmg.lo_freqs,1)),...
        squeeze(mean(cmg.hi_freqs(:,:,1),1)),...
        squeeze(metric)',...
        'linestyle','none');
    caxis([col_levels(1) col_levels(end)]);
    colormap('hot');
    xlabel('Modulating Frequency (Hz)');
    ylabel('Modulated Frequency (Hz)');
    colorbar();
else

    if nargin < 5 || isempty(col_levels)
        % Set colour levels based on the data
        col_levels = 0:max(max(metric))/24:max(max(metric));
    end

    if nargin < 4 || isempty(log_high)
        % Default to not using a log scale for the high frequency
        log_high = false;
    end

    %% Extract the relevant parts of the spectrum
    lo_start = mean(cmg.lo_freqs(:,1));
    lo_stop = mean(cmg.lo_freqs(:,end));
    hi_start = mean(cmg.hi_freqs(:,1,1));
    hi_stop = mean(cmg.hi_freqs(:,end,end));

    lo_freq = freq_vect(lo_start < freq_vect & freq_vect < lo_stop);
    lo_pow = spectrum(lo_start < freq_vect & freq_vect < lo_stop);

    hi_freq = freq_vect(hi_start < freq_vect & freq_vect < hi_stop);
    hi_pow = spectrum(hi_start < freq_vect & freq_vect < hi_stop);

    %% Plot modulating frequency underneath the main comodulogram
    figure;
    axes('position',[.3,.1,.6,.2])
    plot(lo_freq,lo_pow)
    xlim([lo_start lo_stop]);
    xlabel('Modulating Frequency (Hz)');

    %% Plot the modulated frequency spectrum vertically to the left of the
    %  comodulogram
    axes('position',[.1,.3,.2,.6])
    plot(hi_freq,hi_pow)
    ylim([0 max(max(hi_pow))*1.2])
    view(-270, 90);
    set(gca, 'xdir', 'reverse');
    xlim([hi_start hi_stop]);
    xlabel('Modulated Frequency (Hz)');
    % Use log scale if requested
    if log_high == true; set(gca, 'YScale','log'); end

    %% Plot the main comodulogram
    axes('position',[.3,.3,.6,.6])
    contourf(squeeze(mean(cmg.lo_freqs,1)),...
        squeeze(mean(cmg.hi_freqs(:,:,1),1)),...
        squeeze(metric)',...
        'linestyle','none');
    caxis([col_levels(1) col_levels(end)]);
    colormap('hot');
    colorbar('position',[.92,.3,.03,.6]);

end
