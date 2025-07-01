function Rtrn_mat = integrateTrialsAgainstTests(Trials,Tests,dim)

if isa(Trials,"cell")
    if dim == 1
        Rtrn_mat = zeros(length(Tests),size(Trials{1},2),length(Trials));
        for k = 1:length(Tests)
            for j = 1:length(Trials)
                Rtrn_mat(k,:,j) = sum(Tests{k} .* Trials{j},1);
            end
        end
    elseif dim ==2
        Rtrn_mat = zeros(size(Trials{1},1),length(Tests),length(Trials));
        for k = 1:length(Tests)
            for j = 1:length(Trials)
                Rtrn_mat(:,k,j) = sum(Tests{k} .* Trials{j},2);
            end 
        end
    end
elseif length(size(Trials)) == 2
     if dim == 1
        Rtrn_mat = zeros(length(Tests),size(Trials,2));
        for k = 1:length(Tests)
            for j = 1:size(Trials,3)
                Rtrn_mat(k,j) = sum(Tests{k} .* Trials(:,j),1);
            end
        end
    elseif dim ==2
        Rtrn_mat = zeros(size(Trials,1),length(Tests));
        for k = 1:length(Tests)
            for j = 1:size(Trials,3)
                Rtrn_mat(k,j) = sum(Tests{k} .* Trials(:,j),2);
            end 
        end
    elseif dim == 3
       
        %opt for boundary trials
        % Rtrn_mat = zeros(size(length(Tests),size(Trials,2)));
        % for k = 1:length(Tests)
        %     for j = 1:size(Trials,2)
        %         Rtrn_mat(k,j) = sum(Tests{k} .* Trials(:,j)');
        %     end
        % end
        L = 1;
            Rtrn_mat = zeros(length(Tests),L);
            for k = 1:length(Tests)                
                Rtrn_mat(k) = sum(Tests{k} .* Trials);                
            end
    end

else
    if dim == 1
        Rtrn_mat = zeros(length(Tests),size(Trials,2),size(Trials,3));
        for k = 1:length(Tests)
            for j = 1:size(Trials,3)
                Rtrn_mat(k,:,j) = sum(Tests{k} .* Trials(:,:,j),1);
            end
        end
    elseif dim ==2
        Rtrn_mat = zeros(size(Trials,1),length(Tests),size(Trials,3));
        for k = 1:length(Tests)
            for j = 1:size(Trials,3)
                Rtrn_mat(:,k,j) = sum(Tests{k} .* Trials(:,:,j),2);
            end 
        end
    elseif dim == 3
        %opt for boundary trials
        if ndims(Trials) ==3
            
            Trials = squeeze(Trials);
            L = size(Trials,2);
            Rtrn_mat = zeros(size(length(Tests),L));
            for k = 1:length(Tests)
                for j = 1:L
                    Rtrn_mat(k,j) = sum(Tests{k} .* Trials(:,j)');
                end
            end
        else
            
            L = 1;
            Rtrn_mat = zeros(size(length(Tests),L));
            for k = 1:length(Tests)                
                Rtrn_mat(k) = sum(Tests{k} .* Trials);                
            end
        end
        
    end
end