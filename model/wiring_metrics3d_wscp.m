
function [cost, rout_eff, st_cap, R, vol] = wiring_metrics3d_wscp (trees,expscp)
% Calculate R (arbor radius), connectivity repertoire aka storage capacity, routing efficiency and cost
% last two are with the same formula as 2D trees

% R^2 = (1/L^2)sum(i=1, K){sum(j=1 j=!i, K){dli*dlj*(ri-rj)^2}}
% S ~= S0+(L/a)*ln(1-(R/l)) - (l^2)/La - (L^2)/(R^2)

tree = trees;
while(numel(tree)==1)
     tree=tree{1}; 
end
stats = stats_tree(tree,[],[],'-x');
L = stats.gstats.len;
for t = 1:length(tree)
    K = length(len_tree(tree{1, t}));
    dl = len_tree(tree{1, t});
    r_x = tree{1, t}.X;
    r_y = tree{1, t}.Y;
    r_z = tree{1, t}.Z;
    soma = 0;
    for i = 1:K
        for j = 1:K
            if j ~= i
                soma = soma + dl(i)*dl(j)*((r_x(i) - r_x(j))^2 + (r_y(i) - r_y(j))^2 + (r_z(i) - r_z(j))^2);
            elseif i, j = K;
                R2 = (1/L(t)^2)*soma;
                R(t) = sqrt(R2);
            end
        end
    end
end

a = 4;      % [um] Later try to obtain the value by the formula??
s = 2;      % [um] Spine reach length
spd = expscp*2;        % [spines/um]2 Based on WT valued Dierssen 2003
rhoa = 8.15;   % [um/um3] Axon length per unit of volume Figure 4 D https://www.researchgate.net/publication/326139871_The_effects_of_aging_on_neuropil_structure_in_mouse_somatosensory_cortex-A_3D_electron_microscopy_analysis_of_layer_1

tree = trees;
while(numel(tree)==1)
     tree=tree{1}; 
end
%stats = stats_tree(tree, [], [], '-w -s');
%L = stats.gstats.len;
branch_num = stats.gstats.bpoints;
mpeucl = stats.gstats.mpeucl;
for t = 1:length(tree)
    x = tree{1, t}.X;
    y = tree{1, t}.Y;
    z = tree{1, t}.Z;
    P = [x, y, z];
    [k_f, v_f] = boundary(P, 0.7);
    vol(t) = v_f;
    cv(t) = std(stats.dstats.blen{t})/mean(stats.dstats.blen{t});
    l = mean(unique(Pvec_tree(tree{1, t}).*T_tree(tree{1, t})));
%     st_cap(t) = (L(t)/a)*log(1 - (R(t)/l) - (l^2)/(L(t)*a) - (L(t)^2)/R(t)^2);
    %st_cap(t) = (L(t)/a)*log(1 - (R(t)/l)) - (l^2)/(L(t)*a) - expscp*(L(t)^2)/R(t)^2;
    N=s*L(t)*rhoa;
    r=L(t)*spd(t);
    rN = r/N;
    H=(rN*log(1/rN)+(1-rN)*log(1/(1-rN)));
    st_cap(t) =  N*H - s*(L(t)^2)*2/R(t)^2 + (L(t)/a)*(1+log(1 - (R(t)/l))) - (l^2)/(L(t)*a);
    rout_eff(t) = 0;
    cost(t) = stats.gstats.len(t);
end

st_cap = real(st_cap);
end