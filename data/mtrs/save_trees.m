trees = load_tree('wt_aligned.mtr');

for i = 1:length(trees{1})
    swc_tree(trees{1}{i},['WT_2mo_',num2str(i),'.swc']) 
end