function Trials_cell = vectorizeTests(Tests,x)

Trials_cell = cell(length(Tests),1);
dx = mean(diff(x));
for i = 1:length(Tests)
    % Temp_Mat = zeros(size(U_noisy));
    Trials_cell{i}= dx*Tests{i}(x);
end
end