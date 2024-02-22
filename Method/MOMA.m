% load model
nnewmodel=iAQY970

% Obtain a list of reaction IDs for the model
reactionIDs = nnewmodel.rxns;

% Initialize an array to store the values of Fixed ammonia exchange rate and Symbiotic production rate obtained from each loop
f_505_values = zeros(length(reactionIDs), 1);
s_508_values = zeros(length(reactionIDs), 1);

% Initialize an array to store results without records
missingRecords = cell(length(reactionIDs), 1);

for i = 1:length(reactionIDs)
    try
        nnewmodel2 = nnewmodel;

        % Change reaction limits
        nnewmodel2 = changeRxnBounds(nnewmodel2, reactionIDs{i}, 10, 'l');

        % Perform linear MOMA analysis
        [solutionDel, ~, ~, ~] = linearMOMA(nnewmodel1, nnewmodel2, 'max', true);

        % 获取第505个元素
        f_505 = solutionDel.x(505);
        s_508 = solutionDel.x(508);

        % Save the results obtained from each cycle
        f_505_values(i) = f_505;
        s_508_values(i) = s_508;

        % Output the result obtained from each cycle
        fprintf('Reaction: %s, f_505: %f, x_505: %f\n', reactionIDs{i}, f_505, s_508);

        % Check if there are records
        if isnan(f_505) || isnan(s_508)
            missingRecords{i} = reactionIDs{i};
        end
    catch
        % Deal with exceptions (skip current loop)
        fprintf('Error occurred for Reaction: %s. Skipping...\n', reactionIDs{i});
        continue;
    end
end

% Clear cells without records
missingRecords = missingRecords(~cellfun('isempty', missingRecords));

% Display reaction IDs without records
disp('Reaction IDs without records:');
disp(missingRecords);

% 保存 f_505 和 x_505 的值
save('f_s_values.mat', 'f_505_values', 's_508_values');
