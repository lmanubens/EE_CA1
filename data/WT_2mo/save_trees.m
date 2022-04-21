trees = load_tree('wteene_aligned.mtr');

for i = 1:length(trees{1})
    disp(trees{1}{i}.name)
    swc_tree(trees{1}{i},['WT_EE_6mo_',trees{1}{i}.name(1:2),'_',num2str(i),'.swc']) 
end