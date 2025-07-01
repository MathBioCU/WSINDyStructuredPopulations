function [E2,Einfty,TPR,true_w] = computePerformance2(w_pde,w_ode,...
                TransportTrials,SourceTrials,BoundaryTrials,...
                True_pde,True_ode,...
                TransportTags,SourceTags,BoundaryTags,TrueTags)



if ~isempty(TransportTrials) 
    w = [w_pde;w_ode];

    w_trans = w_pde(1:length(TransportTrials));
    w_source = w_pde(length(TransportTrials)+1:end);
    True_trans = True_pde(1:length(TrueTags{1}));
    True_source = True_pde(length(TrueTags{1})+1:end);
    
    true_trans_w = zeros(size(w_trans));
    
    for i = 1:length(TrueTags{1})
    true_trans_w([TransportTags{:}] == TrueTags{1}(i)) = True_trans(i);
    end
    true_source_w = zeros(size(w_source));
    for i = 1:length(TrueTags{2})
    true_source_w([SourceTags{:}] == TrueTags{2}(i)) = ...
        True_source(i);
    end
    true_pde_w = [true_trans_w;true_source_w];

    true_boundary_w = zeros(size(w_ode));
    for i = 1:length(TrueTags{3})
    true_boundary_w([BoundaryTags{:}] == TrueTags{3}(i)) = True_ode(i);
    end
    true_w = [true_pde_w;true_boundary_w];
    
    E2 =norm(w-true_w)/norm(true_w);
    Einfty = max(abs(w(true_w ~= 0)-true_w(true_w ~=0))./abs(true_w(true_w ~=0)));
    TP = sum(and(w,true_w)); 
    FP = sum(and(w,~true_w)); 
    FN = sum(and(~w,true_w));
    TPR = TP/(TP+FP+FN);

    


else %Assumed Age-strucutre case
    
    w_source = w_pde;
    True_source = True_pde;
    
    
    true_source_w = zeros(size(w_source));
    for i = 1:length(TrueTags{1})
    true_source_w([SourceTags{:}] == TrueTags{1}(i)) = True_source(i);
    end
    true_pde_w = [true_source_w];
    
    
    E2_pde =norm(w_pde-true_pde_w)/norm(true_pde_w);
    Einfty_pde = max(abs(w_pde(true_pde_w ~= 0)-true_pde_w(true_pde_w ~=0))./abs(true_pde_w(true_pde_w ~=0)));
    TP = sum(and(w_pde,true_pde_w)); FP = sum(and(w_pde,~true_pde_w)); FN = sum(and(~w_pde,true_pde_w));
    TPR_pde = TP/(TP+FP+FN);
    
    true_boundary_w = zeros(size(w_ode));
    for i = 1:length(TrueTags{2})
    true_boundary_w([BoundaryTags{:}] == TrueTags{2}(i)) = True_ode(i);
    end
    E2_ode = norm(w_ode - true_boundary_w,2)/norm(true_boundary_w,2);
    Einfty_ode = max(abs(w_ode(true_boundary_w ~=0) - true_boundary_w(true_boundary_w~=0))./abs(true_boundary_w(true_boundary_w~=0)));
    TP = sum(and(w_ode,true_boundary_w));  FP = sum(and(w_ode,~true_boundary_w)); FN = sum(and(~w_ode,true_boundary_w));
    TPR_ode = TP/(TP+FP+FN);


end