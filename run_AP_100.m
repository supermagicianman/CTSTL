
cho_set_record = zeros(100, 6); 
model = 'AP';
load_system(model);
closeLoop = 1;
configSet = getActiveConfigSet(model);

for j = 1:100

    cho_set = randi(100, 6, 1); 
    disp(cho_set); 
    cho_set_record(j, :) = cho_set; 
    
    params = get_param(configSet, 'ObjectParameters');

    % Loop through each parameter and apply it to the model
    paramNames = fieldnames(params);
    for i = 1:length(paramNames)
        try
            % Get the value of the parameter
            paramValue = get_param(configSet, paramNames{i});
            
            % Apply the parameter to the model
            set_param(model, paramNames{i}, paramValue);
        catch
            % Skip parameters that cannot be set directly
            warning('Skipping parameter: %s', paramNames{i});
        end
    end

    out_data = sim(model); 
    
    filename = sprintf('./AP_trace/out-%d.mat', j);
    save(filename, 'out_data');
    
    data = load('out.mat');
    logs = data.logsout;
    signal = logs.get(1);
    signal_name = signal.Name;  % Signal Name
    
    % Get time series data
    ts_data = signal.Values;  % Time series object
    
    % Convert to table
    T = table(ts_data.Time, ts_data.Data);
    
    % Define output filename
    filename = sprintf('./AP_trace/trace-%d.csv', j);
    
    % Write to CSV
    writetable(T, filename);
    
end

writecell(cho_set_record, 'cho_set_record.csv')

