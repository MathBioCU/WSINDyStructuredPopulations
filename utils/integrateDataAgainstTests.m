function Rtrn_mat = integrateDataAgainstTests(Data, Testscell,dim)

if dim == 1
    Rtrn_mat = zeros(length(Testscell), size(Data,2));
    for k = 1:length(Testscell)
        Rtrn_mat(k,:) = sum(Testscell{k} .* Data,1);
    end

elseif dim == 2
    Rtrn_mat = zeros(size(Data,1), length(Testscell));
    for k = 1:length(Testscell)
        Rtrn_mat(:,k) = sum(Testscell{k} .* Data,2);
    end
end