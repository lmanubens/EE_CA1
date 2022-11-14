% Let's load first the trees to be simulated
wttrees = load_tree('wt_aligned.mtr');%{sample_tree};   % load sample tree
tgtrees = load_tree('tg_aligned.mtr');
wteetrees = load_tree('wtee_aligned.mtr');%{sample_tree};   % load sample tree
tgeetrees = load_tree('tgee_aligned.mtr');
wtnenetrees = load_tree('wtnene_aligned.mtr');%{sample_tree};   % load sample tree
tgnenetrees = load_tree('tgnene_aligned.mtr');
wteenetrees = load_tree('wteene_aligned.mtr');%{sample_tree};   % load sample tree
tgeenetrees = load_tree('tgeene_aligned.mtr');

% Mature/mushroom spines
% wttrees = wttrees{1};
% wtsyn = 1.3/1.3;
% tgtrees = tgtrees{1};
% tgsyn = 1/1.3;
% wteetrees = wteetrees{1};
% wteesyn = 1.5/1.3;
% tgeetrees = tgeetrees{1};
% tgeesyn = 1.1/1.3;
% wtnenetrees = wtnenetrees{1};
% wtnenesyn = 1.3/1.3;
% tgnenetrees = tgnenetrees{1};
% tgnenesyn = 1/1.3;
% wteenetrees = wteenetrees{1};
% wteenesyn = 1.15/1.3;
% tgeenetrees = tgeenetrees{1};
% tgeenesyn = 0.9/1.3;

% Mushroom spines
% wttrees = wttrees{1};
% wtsyn = .75/.75;
% tgtrees = tgtrees{1};
% tgsyn = .35/.75;
% wteetrees = wteetrees{1};
% wteesyn = .7/.75;
% tgeetrees = tgeetrees{1};
% tgeesyn = .42/.75;
% wtnenetrees = wtnenetrees{1};
% wtnenesyn = .75/.75;
% tgnenetrees = tgnenetrees{1};
% tgnenesyn = .5/.75;
% wteenetrees = wteenetrees{1};
% wteenesyn = .6/.75;
% tgeenetrees = tgeenetrees{1};
% tgeenesyn = .4/.75;

% % Only morphology
wttrees = wttrees{1};
wtsyn = 1;
tgtrees = tgtrees{1};
tgsyn = 1.07;
wteetrees = wteetrees{1};
wteesyn = 1.21;
tgeetrees = tgeetrees{1};
tgeesyn = .93;
wtnenetrees = wtnenetrees{1};
wtnenesyn = 1.07;
tgnenetrees = tgnenetrees{1};
tgnenesyn = 1.07;
wteenetrees = wteenetrees{1};
wteenesyn = 1.18;
tgeenetrees = tgeenetrees{1};
tgeenesyn = .89;

tree = {wttrees{:},tgtrees{:},wteetrees{:},tgeetrees{:},wtnenetrees{:},tgnenetrees{:},wteenetrees{:},tgeenetrees{:}};

expscp = [repmat(wtsyn,length(wttrees),1),
    repmat(tgsyn,length(tgtrees),1),
    repmat(wteesyn,length(wteetrees),1),
    repmat(tgeesyn,length(tgeetrees),1),
    repmat(wtnenesyn,length(wtnenetrees),1),
    repmat(tgnenesyn,length(tgnenetrees),1),
    repmat(wteenesyn,length(wteenetrees),1),
    repmat(tgeenesyn,length(tgeenetrees),1)
    ];

[cost, rout_eff, st_cap, R, vol] = wiring_metrics3d_wscp (tree,expscp);

metrics = [reshape(st_cap, [], 1), reshape(cost, [], 1), reshape(R, [], 1), reshape(vol, [], 1)];
metrics = array2table(metrics, 'VariableNames', {'connectivity_repertoire', 'cost', 'R', 'vol'});
writetable(metrics, 'wiring_metrics.csv', 'WriteVariableNames', true);