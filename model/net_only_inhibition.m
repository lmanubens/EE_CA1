%% *T2N simulations*

% Let's load first the trees to be simulated
wttrees = load_tree('wt_aligned.mtr');%{sample_tree};   % load sample tree
tgtrees = load_tree('tg_aligned.mtr');
% this tree has no soma region yet, so it is perfectly suitable to show you
% how new regions can be defined in TREES toolbox. 
% just make the first two nodes from the root "somatic"
% tree{1}.R(1:2) = 3;     % make the first two nodes a new region
% tree{1}.rnames{3} = 'soma';   % name the region "soma"
% tree{1}.D(1:2) = 10;        % increase diameter of the somatic nodes
i=1;
wttrees = wttrees{1};
tgtrees = tgtrees{1};

tree{1} = wttrees{1};
tree{2} = tgtrees{1};
wtsyn = 1;
synfact{1} = 1;
synfact{2} = 1;

% Let's load first the trees to be simulated
wttrees = load_tree('wt_aligned.mtr');%{sample_tree};   % load sample tree
tgtrees = load_tree('tg_aligned.mtr');
wteetrees = load_tree('wtee_aligned.mtr');%{sample_tree};   % load sample tree
tgeetrees = load_tree('tgee_aligned.mtr');
wtnenetrees = load_tree('wtnene_aligned.mtr');%{sample_tree};   % load sample tree
tgnenetrees = load_tree('tgnene_aligned.mtr');
wteenetrees = load_tree('wteene_aligned.mtr');%{sample_tree};   % load sample tree
tgeenetrees = load_tree('tgeene_aligned.mtr');
% 
% Mature/mushroom spines
% wttrees = wttrees{1};
% wtsyn = 1.3/1.3;
% tgtrees = tgtrees{1};
% tgsyn = 1.05/1.3;
% wteetrees = wteetrees{1};
% wteesyn = 1.51/1.3;
% tgeetrees = tgeetrees{1};
% tgeesyn = 1.08/1.3;
% wtnenetrees = wtnenetrees{1};
% wtnenesyn = 1.2/1.3;
% tgnenetrees = tgnenetrees{1};
% tgnenesyn = 1.04/1.3;
% wteenetrees = wteenetrees{1};
% wteenesyn = 1.02/1.3;
% tgeenetrees = tgeenetrees{1};
% tgeenesyn = 0.87/1.3;
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
% 

% B and C spines
wttrees = wttrees{1};
wtsynB = .5025/1.297;
wtsynC = .7945/1.297;
tgtrees = tgtrees{1};
tgsynB = .725/1.297;
tgsynC = .322/1.297;
wteetrees = wteetrees{1};
wteesynB = .7975/1.297;
wteesynC = .708/1.297;
tgeetrees = tgeetrees{1};
tgeesynB = .619/1.297;
tgeesynC = .4595/1.297;
wtnenetrees = wtnenetrees{1};
wtnenesynB = .78/1.297;
wtnenesynC = .42/1.297;
tgnenetrees = tgnenetrees{1};
tgnenesynB = .725/1.297;
tgnenesynC = .3125/1.297;
wteenetrees = wteenetrees{1};
wteenesynB = .671/1.297;
wteenesynC = .35/1.297;
tgeenetrees = tgeenetrees{1};
tgeenesynB = 0.6125/1.297;
tgeenesynC = 0.255/1.297;

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
% 
% % All spines
% wttrees = wttrees{1};
% wtsyn = 1;
% tgtrees = tgtrees{1};
% tgsyn = 1.07;
% wteetrees = wteetrees{1};
% wteesyn = 1.21;
% tgeetrees = tgeetrees{1};
% tgeesyn = .93;
% wtnenetrees = wtnenetrees{1};
% wtnenesyn = 1.07;
% tgnenetrees = tgnenetrees{1};
% tgnenesyn = 1.06;
% wteenetrees = wteenetrees{1};
% wteenesyn = 1.09;
% tgeenetrees = tgeenetrees{1};
% tgeenesyn = .96;

% % Only morphology
% wttrees = wttrees{1};
% wtsyn = 1;
% tgtrees = tgtrees{1};
% tgsyn = 1;
% wteetrees = wteetrees{1};
% wteesyn = 1;
% tgeetrees = tgeetrees{1};
% tgeesyn = .93;
% wtnenetrees = wtnenetrees{1};
% wtnenesyn = 1;
% tgnenetrees = tgnenetrees{1};
% tgnenesyn = 1;
% wteenetrees = wteenetrees{1};
% wteenesyn = 1;
% tgeenetrees = tgeenetrees{1};
% tgeenesyn = 1;

tree = {wttrees{:},tgtrees{:},wteetrees{:},tgeetrees{:},wtnenetrees{:},tgnenetrees{:},wteenetrees{:},tgeenetrees{:}};

treebak=tree;
tree={tree{1:18},tree{1:18},tree{1:18},tree{1:18},tree{1:18},tree{1:18},tree{1:18},tree{1:17}};


%% Initialize parameters and the neuron structure

%
if ~exist('t2n','file')         % checks if T2N is not already on the Matlab path
    if exist('t2n_Tutorial.mlx','file')  % probably you are in the t2n/Tutorials folder
        addpath(genpath(fileparts(pwd)));
    else
        error('Please run the "t2n_runthisAfterUnzip.m" script in the T2N folder or add the T2N folder including subfolders to the Matlab path!')
    end
end
if ~exist('load_tree','file')
    error('Please run the "start_trees.m" script in the TREES folder or add that folder including subfolders to the Matlab path!')
end


t2n_initModelfolders(pwd);                                  % initialize model folder hierarchy in current folder

neuron = {};                                                % clear neuron structure


for i = 1:length(tree)                                                                   % preconfiguration to loop through several morphologies, now only one morphology is considered
    neuron{i} = [];
%     neuron.mech{t}.all.pas = struct('g',0.0003,'Ra',100,'e',-80,'cm',1);    % add passive channel to all regions and define membrane capacity [�F/cm�], cytoplasmic resistivity [Ohm cm] and e_leak [mV]
%     neuron.mech{t}.soma.pas = struct('g',0.0004,'Ra',100,'e',-80,'cm',1);    % add passive channel to somatic regions (will overwrite the "all" definition) and define membrane capacity [�F/cm�], cytoplasmic resistivity [Ohm cm] and e_leak [mV]
%     neuron.mech{t}.all.k_ion.ek = -90;
%     neuron.mech{t}.all.na_ion.ena = 50;
%     neuron.mech{t}.soma.hh = struct('gnabar',0.25,'gkbar',0.036,'gl',0);     % add Hodgkin-Huxley Sodium channel only to soma

    neuron{i}.experiment = ['test_',num2str(i)];                                                 % give your simulation run a name (necessary when using advanced t2n protocols)

    neuron{i}.params.v_init = -80;                                 % starting membrane potential [mV] of all cells
    neuron{i}.params.dt = 0.025;                                   % integration time step [ms]
    neuron{i}.params.tstop = 1000;%5000;%300;                                  % stop simulation after this (simulation) time [ms]
    neuron{i}.params.prerun = -400;                                % add a pre runtime [ms] to let system settle
    neuron{i}.params.celsius = 10;                                 % temperature [celsius]
    neuron{i}.params.nseg = 'dlambda';                             % the dlambda rule is used to set the amount of segments per section. Alternatively, a fixed number can be entered or a string 'EachX' with X being the interval between segments [micron]
    neuron{i}.params.accuracy = 0;                                 % optional argument if number of segments should be increased in the soma and axon

    neuron{i}.mech{1} = [];
    
    strct = struct;
    
    mechdata = readtable('mechanisms.csv');
    for ii = 1:size(mechdata,1)
        istrct = table2struct(mechdata(ii,:));
        if ~isfield(strct,(istrct.Region)) %&& ~isfield(strct.(istrct.Region),(istrct.Mechanism))
            strct.(istrct.Region).(istrct.Mechanism) = struct(istrct.Parameter,istrct.Value);
        else
            [strct.(istrct.Region).(istrct.Mechanism).(istrct.Parameter)] = istrct.Value;
        end
    end
    neuron{i}.mech{1} =  t2n_catStruct(strct);
    
    neuron{i}.record{1}.cell = struct('node',1,'record','v');                   % record voltage "v" from first node (i.e. the soma)
end
%% Loading morphologies


for i = 1:length(tree)
    tree{i} = resample_tree(tree{i},2);
    tree{i} = quaddiameter_tree(tree{i});
    tree{i} = soma_tree(tree{i});
    tree{i} = rot_tree(tree{i},[],'-m3dY');
    tree{i}.rnames{1} = 'soma';
    tree{i}.rnames{2} = 'adend';
    tree{i}.rnames{3} = 'trunk';
    dists = eucl_tree(tree{i});
    tree{i}.R(find(dists>500))=2;
    tree{i}.R(find(dists<=500))=2;
    tree{i}.R(find(dists<=10))=2;
    plens = PL_tree(tree{i});
    tipid = find(plens==max(plens));
    tree{i}.R(tipid) = 3;
    parents = idpar_tree(tree{i});
    walkid = parents(tipid);
    tree{i}.R(walkid) = 3;
    while parents(walkid)~= 1
        walkid = parents(walkid);
        tree{i}.R(walkid) = 3;
    end
    for j = 1:length(tree{i}.D)
        if (tree{i}.D(j) > 10) tree{i}.R(j) = 1; end
    end
%     figure;
%     plot_tree(tree{i},tree{i}.R);
%     colorbar   % plot tree with region code
%     axis off
    
    widthfac=0.65;
    plen = Pvec_tree(tree{i});                                         % get path lengths to root for each node
    vec1 = double(plen<50*widthfac).*0.1.*0.0001;                                       % create the exponentially decaying range vector for the hh sodium conductance
    vec2 = double(plen>=50*widthfac).*2.*0.0001;                                       % create the exponentially decaying range vector for the hh sodium conductance
    neuron{i}.mech{1}.range.calH = struct('gcalbar',vec1+vec2); 
    
%     neuron.mech{1}.adend = rmfield(neuron.mech{1}.adend,'calH');
%     neuron.mech{1}.trunk = rmfield(neuron.mech{1}.trunk,'calH');
    
    vec = double(plen>=100*widthfac).*4.*0.0003.*plen./350;
    neuron{i}.mech{1}.range.cat = struct('gcatbar',vec); 
    
%    neuron.mech{1}.adend = rmfield(neuron.mech{1}.adend,'cat');
%    neuron.mech{1}.trunk = rmfield(neuron.mech{1}.trunk,'cat');
    
    vec1 = double(plen>50*widthfac & plen<200*widthfac).*2.*0.9.*1.5.*0.03;
    vec2 = double(plen<=50*widthfac | plen>=200*widthfac).*0.5.*0.9.*1.5.*0.03;
    vec3 = double(plen>50*widthfac & plen<200*widthfac).*5.*0.7.*4.5.*0.0001;
    vec4 = double(plen<=50*widthfac | plen>=200*widthfac).*0.5.*0.7.*4.5.*0.0001;
    vec5 = ones(length(plen),1).*0.00075;
    neuron{i}.mech{1}.range.mykca = struct('gkbar',vec1+vec2); 
    neuron{i}.mech{1}.range.kca = struct('gbar',vec3+vec4,'cac',vec5); 
    
%     neuron.mech{1}.adend = rmfield(neuron.mech{1}.adend,'mykca');
%     neuron.mech{1}.trunk = rmfield(neuron.mech{1}.trunk,'mykca');
%     neuron.mech{1}.adend = rmfield(neuron.mech{1}.adend,'kca');
%     neuron.mech{1}.trunk = rmfield(neuron.mech{1}.trunk,'kca');
    
    vec1 = double(plen>500*widthfac).*1.8e-6.*(1.+3.*500./100);
    
    vec2 = double(plen>300*widthfac).*(-89);
    vec3 = double(plen<=300*widthfac & plen>=100*widthfac).*(-81-8*(plen-100)/200);
    vec4 = double(plen<100*widthfac).*(-81);
    
    vec5 = double(plen>300*widthfac).*7.*0.0005.*4;
    vec6 = double(plen<=300*widthfac & plen>=100*widthfac).*7*0.0005.*(1+plen./100);
    
    vec8 = double(plen<=100*widthfac).*7.*0.0005.*(1+plen./100);
    %vec9 = ones(length(plen),1).*-80;
    
%     t2n_getMech(neuron{i},tree{i},'gkabar_kad'); colorbar
    
    neuron{i}.mech{1}.range.h = struct('gbar',vec1,'vhalf',vec2+vec3+vec4); 
    neuron{i}.mech{1}.range.kad = struct('gkabar',vec5+vec6); 
    neuron{i}.mech{1}.range.kap = struct('gkabar',vec8);%,'ek',vec9); 
    
%     neuron.mech{1}.adend = rmfield(neuron.mech{1}.adend,'h');
%     neuron.mech{1}.trunk = rmfield(neuron.mech{1}.trunk,'h');
%     neuron.mech{1}.adend = rmfield(neuron.mech{1}.adend,'kad');
%     neuron.mech{1}.trunk = rmfield(neuron.mech{1}.trunk,'kad');
%     neuron.mech{1}.adend = rmfield(neuron.mech{1}.adend,'kap');
%     neuron.mech{1}.trunk = rmfield(neuron.mech{1}.trunk,'kap');
    
    tree{i} = t2n_writeTrees(tree{i},[],fullfile(pwd,['test_',num2str(i),'.mtr']));                 % transform tree to NEURON morphology (.hoc file). This only has to be done once for each morphology
    
%     t2n_getMech(neuron{i},tree{i},'gkabar_kad'); colorbar

end



%% Small network with Poissonian input Str Lacunosum

data = [];
synIDs = {};
synIDs2 = {};
spiketimes={};
numspikes={};
totspikeswt=zeros(1,1000);%5000);%300);
totspikestg=zeros(1,1000);%5000);%300);
freq=100;
    
% [spikeMat,tvec] = t2n_poissonSpikeGen(freq,neuron{1}.params);                     % generate a 100 Hz Poissonian spike distribution


% tree{37}=tree{1};

for i=1:37%143*4%length(tree)%37*4
%      [spikeMat,tvec] = t2n_poissonSpikeGen(freq,neuron{1}.params);                     % generate a 100 Hz Poissonian spike distribution

    trees = {};
    trees{1} = tree{ceil(i/1)};
%     trees{2} = struct('artificial','VecStim');                             % add a VecStim as an artificial cell. This class can have a vector played to it defining its activity
    trees{2} = struct('artificial','NetStim','start',10,'interval',(1/freq)*1000,'number',ceil(freq*(neuron{1}.params.tstop-20)/100)); % add a NetStim as an artificial cell and define the start (10 ms) the interval (15 ms) and the number (10) of spikings
    trees{3} = struct('artificial','IntFire2','taus',15); % add a IntFire2 as an artificial cell and define its tau decay time
    trees{4} = struct('artificial','IntFire2','taus',15); % add a IntFire2 as an artificial cell and define its tau decay time
    trees = t2n_writeTrees(trees,[],fullfile(pwd,'inh_net_ii.mtr'));                                    % tree morphologies are rewritten because this NetStim might be not written yet

    nneuron = neuron{ceil(i/1)};                                                         % copy standard neuron structure
%     nneuron.params.parallel = 4;                                              % enable parallel Neuron and set the number of used cores to 4
    nneuron.params.cvode = 1;                                                 % enable the variable time step method to fully profit from parallel Neuron (be reminded that out.t is a cell then, as time steps can be different between cells)
    if i/1<=18
%       synfac=wtsyn;
        vglutvgat=1.15;
        synfacB=wtsynB;
        synfacC=wtsynC;
        radlac=240;
    elseif i/1>18 && i/1<=37
%             synfac=tgsyn;
        vglutvgat=.7;
        synfacB=tgsynB;
        synfacC=tgsynC;
        radlac=181;
    elseif i/1>37 && i/1<=56
%             synfac=wteesyn;
        vglutvgat=1.4;
        synfacB=wteesynB;
        synfacC=wteesynC;
        radlac=234;
    elseif i/1>57 && i/1 <=71
%             synfac=tgeesyn;
        vglutvgat=1.15;
        synfacB=tgeesynB;
        synfacC=tgeesynC;
        radlac=220;
    elseif i/1>71 && i/1 <=89
%             synfac=wtnenesyn;
        vglutvgat=1.05;
        synfacB=wtnenesynB;
        synfacC=wtnenesynC;
        radlac=248;
    elseif i/1>89 && i/1 <=107
%             synfac=tgnenesyn;
        vglutvgat=.6;
        synfacB=tgnenesynB;
        synfacC=tgnenesynC;
        radlac=207;
    elseif i/1>107 && i/1 <=125
%             synfac=wteenesyn;
        vglutvgat=1.15;
        synfacB=wteenesynB;
        synfacC=wteenesynC;
        radlac=241;
    elseif i/1>125 && i/1 <=143
%             synfac=tgeenesyn;
        vglutvgat=.7;
        synfacB=tgeenesynB;
        synfacC=tgeenesynC;
        radlac=206;
    end
    radlac=0.65*radlac;
%         radlac=130;
    %inds = find(50<euclen<150);
%     inds = find(30<trees{1}.Y<radlac);
    inds = find(30<trees{1}.Y);
    lens = len_tree(trees{1});
    lenfac = sum(lens(inds))/1.0611e+03;
    LOs = LO_tree(trees{1});
    probvec = 0.6 + 0.4*(LOs(inds)/max(LOs(inds)));
%         synIDs{i} = randsample(inds,round(sum(lens(inds))*2*0.15*synfac),true,probvec); 
    synIDs{i} = randsample(inds,round(sum(lens(inds))*2*0.075*synfacB),true,probvec);
    synIDs2{i} = randsample(inds,round(sum(lens(inds))*2*0.075*synfacC),true,probvec);
    
    tags=cell(length(synIDs{i})+length(synIDs2{i}),1);
    for j = 1:length(synIDs{i} )+length(synIDs2{i})
        tags{j} = ['EcxSyn_',num2str(j)];
    end
    tags=string(tags);
    nneuron.pp{1}.Exp2Syn(1) = struct('node',[synIDs{i};synIDs2{i}],'tau1',0.2,'tau2',2,'e',0,'tag',tags);% add an Exp2Syn at this location with 0.2 ms rise and 2.5 ms decay time
    
%     tags=cell(length(synIDsExc),1);
%     for j = 1:length(synIDsExc)
%         tags{j} = ['EcxSyn_',num2str(j)];
%     end
%     tags=string(tags);
%     nneuron.pp{1}.Exp2Syn(1) = struct('node',synIDsExc,'tau1',0.2,'tau2',2,'e',0,'tag',tags);% add an Exp2Syn at this location with 0.2 ms rise and 2.5 ms decay time
    synIDsInh = 1;                                                            % choose root as point of inhibitory synapse
    nneuron.pp{1}.Exp2Syn(2) = struct('node',synIDsInh,'tau1',0.5,'tau2',5,'e',-80,'tag',string({'InhSyn'}));% add an inhibitory (e=-80) Exp2Syn at this location with 0.2 ms rise and 2.5 ms decay time
    %nneuron.record{1}.cell.node = cat(1,1,synIDsExc);                         % record somatic voltage and voltage at synapse
    nneuron.record{1}.cell.node = [1];                            % record somatic voltage and voltage at synapse
%     nneuron.record{3}.cell.node = [1];                            % record somatic voltage and voltage at synapse
%     nneuron.record{1}.Exp2Syn(1) = struct('node',synIDs{i},'record','i');        % record synaptic current
%     nneuron.record{1}.Exp2Syn(2) = struct('node',synIDsInh,'record','i');        % record synaptic current
    nneuron.con(1) = struct('source',struct('cell',2),'target',struct('cell',1,'pp','Exp2Syn','node',synIDs{i}),'delay',.1,'threshold',0.5,'weight',.00018);  % connect the NetStim (cell 2) with the target (point process Exp2Syn of cell 1 at node specified in synIDs), and add threshold/weight and delay of the connection (NetStim parameters)
    nneuron.con(2) = struct('source',struct('cell',2),'target',struct('cell',1,'pp','Exp2Syn','node',synIDs2{i}),'delay',.1,'threshold',0.5,'weight',.00012);  % connect the NetStim (cell 2) with the target (point process Exp2Syn of cell 1 at node specified in synIDs), and add threshold/weight and delay of the connection (NetStim parameters)

%     nneuron.con(1) = struct('source',struct('cell',2),'target',struct('cell',1,'pp','Exp2Syn','node',synIDsExc),'delay',0.1,'threshold',0.5,'weight',synfac);  % connect the NetStim (cell 2) with the target (point process Exp2Syn of cell 1 at node specified in synIDsExc), and add threshold/weight and delay of the connection (NetStim parameters)
    %nneuron.con(2) = struct('source',struct('cell',2),'target',struct('cell',3),'delay',0.1,'threshold',0.5,'weight',50*lenfac*vglutvgat);  % connect the NetStim (cell 2) with the inhibitory cell to create a feed-forward inhibitory loop, and add threshold/weight and delay of the connection (NetStim parameters)
    synfac = synfacB+synfacC;
%     nneuron.con(3) = struct('source',struct('cell',1,'node',1),'target',struct('cell',3),'delay',0.1,'threshold',0.5,'weight',round(sum(lens(inds))*2*0.075*synfac*.00018/2));  % connect the real cell with the inhibitory cell to create a feed-back inhibitory loop, and add threshold/weight and delay of the connection (NetStim parameters)
    %https://www.sciencedirect.com/science/article/abs/pii/S0306452200004966?via%3Dihub
    %proportions of excitatory and inhibitory synapses
    nneuron.con(3) = struct('source',struct('cell',2),'target',struct('cell',3),'delay',0.1,'threshold',0.5,'weight',round(sum(lens(inds))*2*0.075*synfac*.00018*2));  % connect the real cell with the inhibitory cell to create a feed-back inhibitory loop, and add threshold/weight and delay of the connection (NetStim parameters)
    nneuron.con(4) = struct('source',struct('cell',3),'target',struct('cell',1,'pp','Exp2Syn','node',synIDsInh),'delay',3,'threshold',0.5,'weight',round(sum(lens(inds))*2*0.075*synfac*.0008/15)/1);  % connect the inhibitory cell with the target (point process Exp2Syn of cell 1 at node specified in synIDsInh), and add threshold/weight and delay of the connection (NetStim parameters)
    nneuron.con(5) = struct('source',struct('cell',2),'target',struct('cell',4),'delay',0.1,'threshold',0.5,'weight',round(sum(lens(inds))*2*0.075*synfac*.00018*2));  % connect the real cell with a second inhibitory cell to create a feed-forward inhibitory loop, and add threshold/weight and delay of the connection (NetStim parameters)
    nneuron.con(6) = struct('source',struct('cell',4),'target',struct('cell',3,'pp','Exp2Syn','node',synIDsInh),'delay',3,'threshold',0.5,'weight',round(sum(lens(inds))*2*0.075*synfac*.0008/30)/vglutvgat);  % connect the inhibitory cell with the target inhibitory cell (point process Exp2Syn of cell 1 at node specified in synIDsInh), and add threshold/weight and delay of the connection (NetStim parameters)
    
%     nneuron.con(2) = struct('source',struct('cell',1,'node',1),'target',struct('cell',3),'delay',0,'threshold',0.5,'weight',round(50*synfac*lenfac));  % connect the real cell with the inhibitory cell to create a feed-back inhibitory loop, and add threshold/weight and delay of the connection (NetStim parameters)
%     nneuron.con(3) = struct('source',struct('cell',3),'target',struct('cell',1,'pp','Exp2Syn','node',synIDsInh),'delay',3,'threshold',0.5,'weight',round(50*synfac*lenfac)/vglutvgat);  % connect the inhibitory cell with the target (point process Exp2Syn of cell 1 at node specified in synIDsInh), and add threshold/weight and delay of the connection (NetStim parameters)

%     Needed for Poisson stimulation
%     nneuron.play{2}.cell =  struct('node',1,'play','spike','times',tvec(spikeMat));   % play this spike vector to the VecStim

    out = t2n(nneuron,trees,'-q',['exch_',num2str(i)]);                                    % execute t2n

    neuronid{i} = i;
    [~,ind] = findpeaks(out.record{1}.cell.v{1},'MinPeakHeight',0);
    spiketimes{i} = out.t(ind);
    numspikes{i} = zeros(1,1000);%5000);%300);
    for j=1:length(spiketimes{i})
        numspikes{i}(round(spiketimes{i}(j)))=1;
    end
    if i/1<=18
        totspikeswt=totspikeswt+numspikes{i};
    elseif i/1<=37
        totspikestg=totspikestg+numspikes{i};
    end
%     freqoutp{i} = numel(ind)/0.3;
   
    % plot the result (Vmem at soma and synapse and synaptic current)
%     fig = figure;
%     subplot(4,1,2)
%     hold all
%     plot(out.t{1},out.record{1}.cell.v{1})       % plot time vs voltage at soma
%     plot(out.t{1},out.record{1}.cell.v{synIDsExc})  % plot time vs voltage at dendrite end
%     legend('Loc soma','Loc exc. synapse','Location','northoutside')
%     ylabel('Membrane potential [mV]')
%     subplot(4,1,3)% 
%     % plot the result (Vmem at soma and synapse and synaptic current)
%     fig = figure;
%     subplot(4,1,2)
%     hold all
%     plot(out.t{1},out.record{1}.cell.v{1})       % plot time vs voltage at soma
%     plot(out.t{1},out.record{1}.cell.v{synIDsExc})  % plot time vs voltage at dendrite end
%     legend('Loc soma','Loc exc. synapse','Location','northoutside')
%     ylabel('Membrane potential [mV]')
%     subplot(4,1,3)
%     plot(out.t{1},out.record{1}.Exp2Syn.i{synIDsExc})  % plot time vs synaptic current
%     ylabel('Exc. syn. current [nA]')
%     subplot(4,1,4)
%     plot(out.t{1},out.record{1}.Exp2Syn.i{synIDsInh})  % plot time vs synaptic current
%     ylabel('Inh. syn. current [nA]')
% 
%     subplot(4,1,1)
%     t2n_plotRaster(spikeMat,tvec)
%     title('Spikes NetStim (Poisson distr.)')
%     fig.Position(4) = fig.Position(4) * 2;
%     linkaxes(get(fig,'Children'),'x')
% 
%     plot(out.t{1},out.record{1}.Exp2Syn.i{synIDsExc})  % plot time vs synaptic current
%     ylabel('Exc. syn. current [nA]')
%     subplot(4,1,4)
%     plot(out.t{1},out.record{1}.Exp2Syn.i{synIDsInh})  % plot time vs synaptic current
%     ylabel('Inh. syn. current [nA]')
% 
%     subplot(4,1,1)
%     fig.Position(4) = fig.Position(4) * 2;
%     linkaxes(get(fig,'Children'),'x')
end

csvwrite('numspikes_onlyI.csv',horzcat(numspikes{:})')

%% Plotting
figure
subplot(4,1,[1,3])
t2n_plotRaster(spiketimes)
subplot(4,1,4)
hold on
t2n_plotRaster(10:(1/freq)*1000:1000)%300)
% t2n_plotRaster(spikeMat,tvec)
avact = movmean(totspikeswt,6);
plot(avact)
avactwt = avact;
avact = movmean(totspikestg,6);
plot(avact)
avacttg = avact;
% t2n_plotRaster(spikeMat,tvec)
title('Spikes NetStim (Poisson distr.)')


figure

y = fft(avactwt);
fs = 1000; % sample frequency (Hz)
n = length(avactwt);          % number of samples
f = (0:n-1)*(fs/n);     % frequency range
power = abs(y).^2/n;    % power of the DFT
% powerplt = movmean(power,20);
powerplt=power;
powerpltwt=power;

hold on
plot(f,powerplt)
xlim([1 300])
ylim([0 max(powerplt(50:end))*1.2])

y = fft(avacttg);
fs = 1000; % sample frequency (Hz)
n = length(avacttg);          % number of samples
f = (0:n-1)*(fs/n);     % frequency range
power = abs(y).^2/n;    % power of the DFT
% powerplt = movmean(power,20);
powerplt=power;
xlabel('Frequency')
ylabel('Power')
xlim([1 300])
% ylim([0 max(powerplt(50:end))])
plot(f,powerplt)


cHeader = {'Frequency' 'Power_WT' 'Power_TG'}; %dummy header
commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); %cHeader in text with commas
%write header to file
fid = fopen('pspectrum_onlyI.csv','w'); 
fprintf(fid,'%s\n',textHeader)
fclose(fid)
%write data to end of file
dlmwrite('pspectrum_onlyI.csv',[f',powerpltwt',powerplt'],'-append');
