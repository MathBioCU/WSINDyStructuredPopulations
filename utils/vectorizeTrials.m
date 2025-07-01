function Trials_cell = vectorizeTrials(Trials,t,x,U_noisy,N)

Trials_cell = cell(length(Trials),1);

for i = 1:length(Trials)
    % Temp_Mat = zeros(size(U_noisy));
    Trials_cell{i}= Trials{i}(x,U_noisy,N);
end
end