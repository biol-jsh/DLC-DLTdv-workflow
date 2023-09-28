function lme_output = formatForLME(parameter,speeds,treatments,individual, parameterName)

parameter_pre_table = [];
speed_pre_table = [];
treatment_pre_table = [];
individual_pre_table = [];


for speed_i = 1:size(parameter,1)
    for treatment_i = 1:size(parameter,2)
        
        temp = cell2mat(parameter(speed_i,treatment_i));
        parameter_pre_table = [parameter_pre_table temp];

        temp = (individual(speed_i,treatment_i));
        temp = temp{1};

        individual_pre_table = [individual_pre_table temp];
        speed_pre_table = [speed_pre_table ones(size(temp))*speeds(speed_i)];
        treatment_pre_table = [treatment_pre_table repmat(string(treatments{treatment_i}), 1,numel(temp))];
    end
end

parameter_table = table(parameter_pre_table', speed_pre_table', treatment_pre_table', individual_pre_table', VariableNames=[parameterName,"speed","treatment", "individual"]);
lme_output = fitlme(parameter_table, strcat(parameterName, " ~ treatment + speed + (1|individual)"));
