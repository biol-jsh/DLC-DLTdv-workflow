% piecewise aggregated approximation
function [x_mean, x_std, bin_N, bins, norm_time, x_iqr] = paa(x_input, x_time, output_size, overlap)

bins = cell(1,output_size);
norm_time = nan(1,output_size);
x_mean = nan(1,output_size);
x_std = nan(1,output_size);
x_iqr = nan(1,output_size);
bin_N = nan(1,output_size);

for K = 1:output_size
    bins{K}.norm_time_start = (K-1)/output_size;
    bins{K}.norm_time_mean = (K-.5)/output_size;
    bins{K}.norm_time_end = (K)*1/output_size;
    bins{K}.values = [];
    
    norm_time(K) = bins{K}.norm_time_mean;
    for h=1:size(x_input,1)
        input_h = x_input(h,:);
        time_h = x_time(h,:);
        input_h(isnan(time_h)) = [];
        
        norm_time_h = time_h/max(time_h);
        if(K == 1)
            bins{K}.values = [bins{K}.values ...
                input_h(bins{K}.norm_time_start <= norm_time_h & norm_time_h <= bins{K}.norm_time_end)];
        else
            bins{K}.values = [bins{K}.values ...
                input_h(bins{K}.norm_time_start < norm_time_h & norm_time_h <= bins{K}.norm_time_end)];
        end
    end
    bins{K}.N = numel(bins{K}.values);
    bins{K}.average = mean(bins{K}.values);
    bins{K}.std = std(bins{K}.values);
    bins{K}.IQR = iqr(bins{K}.values);
    
    bin_N(K) = bins{K}.N;
    x_mean(K) = bins{K}.average;
    x_std(K) = bins{K}.std;
    x_iqr(K) = bins{K}.IQR;
end