%% *T2N simulations*

% Note: T2N requires NEURON and TREES to be installed! See Documentation for further information.
% Note for Mac users: Matlab has to be run from a Terminal for T2N to work 
% properly.

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

% Only morphology
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


%% Initialize parameters and the neuron structure
% First we need to be sure that Matlab finds all necessary files located in 
% the T2N and TREES folder. 
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
% 
% Now we need to initialize the t2n parameter structure which includes general 
% NEURON settings. Most of these settings are set to a default value if not explicitly 
% set, but you might want to control most of them

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
    neuron{i}.params.tstop = 300;                                  % stop simulation after this (simulation) time [ms]
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
    tree{i} = quaddiameter_tree(tree{i});
    tree{i} = soma_tree(tree{i});
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
    
    plen = Pvec_tree(tree{i});                                         % get path lengths to root for each node
    vec1 = double(plen<50).*0.1.*0.0001;                                       % create the exponentially decaying range vector for the hh sodium conductance
    vec2 = double(plen>=50).*2.*0.0001;                                       % create the exponentially decaying range vector for the hh sodium conductance
    neuron{i}.mech{1}.range.calH = struct('gcalbar',vec1+vec2); 
    
%     neuron.mech{1}.adend = rmfield(neuron.mech{1}.adend,'calH');
%     neuron.mech{1}.trunk = rmfield(neuron.mech{1}.trunk,'calH');
    
    vec = double(plen>=100).*4.*0.0003.*plen./350;
    neuron{i}.mech{1}.range.cat = struct('gcatbar',vec); 
    
%    neuron.mech{1}.adend = rmfield(neuron.mech{1}.adend,'cat');
%    neuron.mech{1}.trunk = rmfield(neuron.mech{1}.trunk,'cat');
    
    vec1 = double(plen>50 & plen<200).*2.*0.9.*1.5.*0.03;
    vec2 = double(plen<=50 | plen>=200).*0.5.*0.9.*1.5.*0.03;
    vec3 = double(plen>50 & plen<200).*5.*0.7.*4.5.*0.0001;
    vec4 = double(plen<=50 | plen>=200).*0.5.*0.7.*4.5.*0.0001;
    vec5 = ones(length(plen),1).*0.00075;
    neuron{i}.mech{1}.range.mykca = struct('gkbar',vec1+vec2); 
    neuron{i}.mech{1}.range.kca = struct('gbar',vec3+vec4,'cac',vec5); 
    
%     neuron.mech{1}.adend = rmfield(neuron.mech{1}.adend,'mykca');
%     neuron.mech{1}.trunk = rmfield(neuron.mech{1}.trunk,'mykca');
%     neuron.mech{1}.adend = rmfield(neuron.mech{1}.adend,'kca');
%     neuron.mech{1}.trunk = rmfield(neuron.mech{1}.trunk,'kca');
    
    vec1 = double(plen>500).*1.8e-6.*(1.+3.*500./100);
    
    vec2 = double(plen>300).*(-89);
    vec3 = double(plen<=300 & plen>=100).*(-81-8*(plen-100)/200);
    vec4 = double(plen<100).*(-81);
    
    vec5 = double(plen>300).*7.*0.0005.*4;
    vec6 = double(plen<=300 & plen>=100).*7*0.0005.*(1+plen./100);
    
    vec8 = double(plen<=100).*7.*0.0005.*(1+plen./100);
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




%% Synaptic stimulation, Exp2Syn synapse and a NetStim

vsoma = {};
time = {};
data = [];

for freq=10:10:100
    outs = {};
    synIDs = {};
    for i=108:length(tree)
        trees = {};
        trees{1} = tree{i};
        trees{2} = struct('artificial','NetStim','start',10,'interval',(1/freq)*1000,'number',ceil(freq*(neuron{1}.params.tstop-20)/100)); % add a NetStim as an artificial cell and define the start (10 ms) the interval (15 ms) and the number (10) of spikings
        trees = t2n_writeTrees(trees,[],fullfile(pwd,['test',num2str(i),'.mtr']));                                    % tree morphologies are rewritten because this NetStim might be not written yet

        nneuron = neuron{i};                                                         % copy standard neuron structure
        %plen = Pvec_tree(trees{1});                                                % get path length to soma at each node of morphology
        euclen = eucl_tree(trees{1});                                                % get path length to soma at each node of morphology
        %[~,synIDs{i}] = max(plen);                                                   % search for the most far away point in the morpholgy    
        %synIDs{i} = ceil(synIDs{i}/2);
        if i<=18
            synfac=wtsyn;
        elseif i>18 && i<=37
            synfac=tgsyn;
        elseif i>37 && i<=56
            synfac=wteesyn;
        elseif i>57 && i <=71
            synfac=tgeesyn;
        elseif i>71 && i <=89
            synfac=wtnenesyn;
        elseif i>89 && i <=107
            synfac=tgnenesyn;
        elseif i>107 && i <=125
            synfac=wteenesyn;
        elseif i>125 && i <=143
            synfac=tgeenesyn; 
        end
        synIDs{i} = randsample(find(50<euclen<150),round(20*synfac));
        nneuron.pp{1}.Exp2Syn = struct('node',synIDs{i},'tau1',0.2,'tau2',2.5,'e',0);% add an Exp2Syn at this location with 0.2 ms rise and 2.5 ms decay time
        nneuron.record{1}.cell.node = [1,1:10:numel(tree{i}.R),synIDs{i}'];                            % record somatic voltage and voltage at synapse
        nneuron.record{1}.Exp2Syn = struct('node',synIDs{i},'record','i');           % record synaptic current
        %introduce synaptic voltages
        nneuron.con(1) = struct('source',struct('cell',2),'target',struct('cell',1,'pp','Exp2Syn','node',synIDs{i}),'delay',0,'threshold',0.5,'weight',.00375);  % connect the NetStim (cell 2) with the target (point process Exp2Syn of cell 1 at node specified in synIDs), and add threshold/weight and delay of the connection (NetStim parameters)

        out = t2n(nneuron,trees,'-q');                                    % execute t2n
        outs{i} = out;
        
        neuronid{i} = i;
        [~,ind] = findpeaks(out.record{1}.cell.v{1},'MinPeakHeight',0);
        freqout{i} = numel(ind)/0.3;
        freqin{i} = freq;
    end
    data = [data; [cell2mat(neuronid);cell2mat(freqin);cell2mat(freqout)]'];


    % plot the result (Vmem at soma and synapse and synaptic current)
%     figure;
%     subplot(2,1,1)
%     hold all
% %     for i=1:length(outs)
% %         plot(outs{i}.t,outs{i}.record{1}.cell.v{1})       % plot time vs voltage at soma
% %         plot(outs{i}.t,outs{i}.record{1}.cell.v{synIDs{i}(1)})  % plot time vs voltage at dendrite end
% %     end
%     legend('Soma - WT','Synapse - WT','Soma - TG','Synapse - TG')
%     ylabel('Membrane potential [mV]')
%     subplot(2,1,2)
% %     for i=1:length(outs)
% %         plot(outs{i}.t,outs{i}.record{1}.Exp2Syn.i{synIDs{i}(1)})  % plot time vs synaptic current
% %     end
%     ylabel('Synaptic current [nA]')
%     xlabel('Time [ms]')

%     for i = 1:length(tree)
%         plot_syninp(tree{i}, outs{i}, 1, 10/neuron.params.dt);
%     end
end

cHeader = {'neuronid' 'FreqIn' 'FreqOut'}; %dummy header
commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); %cHeader in text with commas
%write header to file
fid = fopen('freqinout.csv','w'); 
fprintf(fid,'%s\n',textHeader)
fclose(fid)
%write data to end of file
dlmwrite('freqinout.csv',data,'-append');

