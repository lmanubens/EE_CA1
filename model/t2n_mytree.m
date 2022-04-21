%% *T2N Tutorial*
% *This script is to show you some applications of t2n. Tutorial A and B initialize 
% the morphology and standard parameters (and need thus be run at least once before 
% starting one of the other tutorials), whereas the other sections define the 
% different simulations.*
% 
% Note: T2N requires NEURON and TREES to be installed! See Documentation for further information.
% Note for Mac users: Matlab has to be run from a Terminal for T2N to work 
% properly.

% Let's load first the trees to be simulated
wttrees = load_tree('wt_aligned.mtr');%{sample_tree};   % load sample tree
tgtrees = load_tree('tg_aligned.mtr');
wteetrees = load_tree('wtee_aligned.mtr');%{sample_tree};   % load sample tree
tgeetrees = load_tree('tgee_aligned.mtr');
wtnenetrees = load_tree('wtnene_aligned.mtr');%{sample_tree};   % load sample tree
tgnenetrees = load_tree('tgnene_aligned.mtr');
wteenetrees = load_tree('wteene_aligned.mtr');%{sample_tree};   % load sample tree
tgeenetrees = load_tree('tgeene_aligned.mtr');

% this tree has no soma region yet, so it is perfectly suitable to show you
% how new regions can be defined in TREES toolbox. 
% just make the first two nodes from the root "somatic"
% tree{1}.R(1:2) = 3;     % make the first two nodes a new region
% tree{1}.rnames{3} = 'soma';   % name the region "soma"
% tree{1}.D(1:2) = 10;        % increase diameter of the somatic nodes
i=1;
wttrees = wttrees{1};
tgtrees = tgtrees{1};
wteetrees = wteetrees{1};
tgeetrees = tgeetrees{1};
wtnenetrees = wtnenetrees{1};
tgnenetrees = tgnenetrees{1};
wteenetrees = wteenetrees{1};
tgeenetrees = tgeenetrees{1};

%tree = {wttrees{:},tgtrees{:},wteetrees{:},tgeetrees{:},wtnenetrees{:},tgnenetrees{:},wteenetrees{:},tgeenetrees{:}};
 tree{1} = wttrees{1};
 tree{2} = tgtrees{1};


%% Tutorial A - initialize parameters and the neuron structure
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

neuron = [];                                                % clear neuron structure
neuron.params.v_init = -80;                                 % starting membrane potential [mV] of all cells
neuron.params.dt = 0.025;                                   % integration time step [ms]
neuron.params.tstop = 300;                                  % stop simulation after this (simulation) time [ms]
neuron.params.prerun = -400;                                % add a pre runtime [ms] to let system settle
neuron.params.celsius = 10;                                 % temperature [celsius]
neuron.params.nseg = 'dlambda';                             % the dlambda rule is used to set the amount of segments per section. Alternatively, a fixed number can be entered or a string 'EachX' with X being the interval between segments [micron]
neuron.params.accuracy = 0;                                 % optional argument if number of segments should be increased in the soma and axon

%
% Now we set the t2n neuron structure which defines all mechanisms, point 
% processes, synapses, connections recordings, protocols etc. that should be used 
% in the simulation. The mechanisms and point processes can be any that is defined 
% in the standard NEURON environment or that has been  written as .mod file and 
% saved to the 'lib_mech' folder in the model folder from which it can be compiled 
% by t2n. 'neuron' is a Matlab structure  if a single simulation should be run.


neuron.experiment = 'test';                                                 % give your simulation run a name (necessary when using advanced t2n protocols)
for t = 1%:length(tree)                                                                   % preconfiguration to loop through several morphologies, now only one morphology is considered
%     neuron.mech{t}.all.pas = struct('g',0.0003,'Ra',100,'e',-80,'cm',1);    % add passive channel to all regions and define membrane capacity [�F/cm�], cytoplasmic resistivity [Ohm cm] and e_leak [mV]
%     neuron.mech{t}.soma.pas = struct('g',0.0004,'Ra',100,'e',-80,'cm',1);    % add passive channel to somatic regions (will overwrite the "all" definition) and define membrane capacity [�F/cm�], cytoplasmic resistivity [Ohm cm] and e_leak [mV]
%     neuron.mech{t}.all.k_ion.ek = -90;
%     neuron.mech{t}.all.na_ion.ena = 50;
%     neuron.mech{t}.soma.hh = struct('gnabar',0.25,'gkbar',0.036,'gl',0);     % add Hodgkin-Huxley Sodium channel only to soma
    
    neuron.mech{t} = [];
    
    strct = struct;
    
    mechdata = readtable('mechanisms.csv');
    for i = 1:size(mechdata,1)
        istrct = table2struct(mechdata(i,:));
        if ~isfield(strct,(istrct.Region)) %&& ~isfield(strct.(istrct.Region),(istrct.Mechanism))
            strct.(istrct.Region).(istrct.Mechanism) = struct(istrct.Parameter,istrct.Value);
        else
            [strct.(istrct.Region).(istrct.Mechanism).(istrct.Parameter)] = istrct.Value;
        end
    end
    neuron.mech{t} =  t2n_catStruct(strct);
    
    neuron.record{t}.cell = struct('node',1,'record','v');                   % record voltage "v" from first node (i.e. the soma)
end
%% Tutorial B - loading a morphology
% First we need a morphology (tree), on which we can add the channels etc. We 
% are using here a sample_tree from the TREES toolbox which loads a small tree. 
% You can also load a different tree by using "tree = load_tree;", or load another 
% example, e.g. "tree = {hsn_tree};"
%

figure;
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
    plot_tree(tree{i},tree{i}.R);
    
%     plen = Pvec_tree(tree{i});                                         % get path lengths to root for each node
%     vec1 = double(plen<50).*0.1.*0.0001;                                       % create the exponentially decaying range vector for the hh sodium conductance
%     vec2 = double(plen>=50).*2.*0.0001;                                       % create the exponentially decaying range vector for the hh sodium conductance
%     neuron.mech{i}.range.calH = struct('gcabar_calH',vec1+vec2); 
%     
%     vec = double(plen>=100).*4.*0.0003.*plen./350;
%     neuron.mech{i}.range.calH = struct('gcatbar_cat',vec); 
%     
%     vec1 = double(plen>50 & plen<200).*2.*0.9.*1.5.*0.03;
%     vec2 = double(plen<=50 | plen>=200).*0.5.*0.9.*1.5.*0.03;
%     vec3 = double(plen>50 & plen<200).*5.*0.7.*4.5.*0.0001;
%     vec4 = double(plen<=50 | plen>=200).*0.5.*0.7.*4.5.*0.0001;
%     neuron.mech{i}.range.mykca = struct('gbar_mykca',vec1+vec2); 
%     neuron.mech{i}.range.kca = struct('gbar_kca',vec3+vec4); 
%     
%     vec1 = double(plen>500).*1.8e-6.*(1.+3.*500./100);
%     
%     vec2 = double(plen>300).*(-89);
%     vec3 = double(plen<=300 & plen>=100).*(-81-8*(plen-100)/200);
%     vec4 = double(plen<100).*(-81);
%     
%     
%     
%     neuron.mech{i}.range.h = struct('gbar_h',vec1); 
%     neuron.mech{i}.range.h = struct('vhalf_h',vec2+vec3+vec4); 
%     neuron.mech{i}.range.kad = struct('gkabar_kad',vec5+vec6+vec7); 
%     neuron.mech{i}.range.kap = struct('gkabar_kap',vec8+vec9+vec10); 
    
end


colorbar   % plot tree with region code
axis off

%tree = tree(1);                                                              % for simplicity, only one tree is looked at (if several exist at all)
tree = t2n_writeTrees(tree,[],fullfile(pwd,'test.mtr'));                 % transform tree to NEURON morphology (.hoc file). This only has to be done once for each morphology
%% Tutorial C - simulation protocol: somatic current injection
% Now we want to do a simple somatic current injection, execute neuron and plot 
% the result. By copying the 'neuron' structure into a new variable 'nneuron' and modify 
% only this one in each protocol, we do not have to rerun the upper sections each time.
%
% nneuron = neuron;                                                           % copy standard neuron structure
% 
% % now we define a 100 ms current injection of 0.6 nA that is started at 
% % time point 50 ms
% nneuron.pp{1}.IClamp = struct('node',1,'times',[50 100 150 200],'amp',[3 0 3 0]);     % add a current clamp electrode to the first node and define stimulation times [ms] and amplitudes [nA]
% nneuron.record{1}.IClamp = struct('node',1,'record','i');                   % record the current of the IClamp just for visualization
% 
% outs={};
% for i=1:length(tree)
%     out = t2n(nneuron,tree{i},'-w-q');                                      % execute t2n and receive output
%     outs{i}=out;
% end
% % 
% % After execution of t2n we can plot the results. T2N returns all recordings 
% % that had previously been defined. We can access them in out.record{1} where 
% % {1} accesses the recordings of the first tree/morphology. Subsequently all recordings 
% % of distributed mechanisms or other cell-specific parameters (e.g. voltage) can 
% % be found in the field  'cell', whereas recorded values from point processes 
% % can be found under the name of that process (e.g. 'IClamp'). Following that 
% % comes the name of the recorded variable and the index to the node at which the 
% % variable has been recorded or the corresponding point process had been placed  
% % (e.g. v{1} for the voltage at the first node). The time vector can be found 
% % in out.t and fits to all recordings, except if the options params.local_dt has 
% % been used which allows NEURON to use different time steps for different trees. 
% % In that case out.t is a cell array of vectors.
% for i = 1:length(tree)
%     plot_syninp(tree{i}, outs{i}, 1, 45/neuron.params.dt);
% end
% 
% figure; 
% subplot(2,1,1)                                                              % make a subplot in the figure
% hold on
% for i=1:length(outs)
%     plot(outs{i}.t,outs{i}.record{1}.cell.v{1})                                         % plot recorded somatic voltage (time vs voltage)
% end
% hold off
% ylim([-90,60])
% xlim([0,nneuron.params.tstop])
% ylabel('Membrane potential [mV]')
% xlabel('Time [ms]')
% subplot(2,1,2)                                                              % make another subplot in the figure
% hold on
% for i=1:length(outs)
%     plot(outs{i}.t,outs{i}.record{1}.IClamp.i{1})                                       % plot electrode current (time vs current)
% end
% hold off
% ylim([0,4])
% xlim([0,nneuron.params.tstop])
% ylabel('Injected current [nA]')
% %% Tutorial D - simulation protocol: several simulations with different injected current amplitudes (f-I relationship)
% % Now we learn something about the real strength of t2n: It can execute simulations 
% % in parallel! All we have have to do is to create a cell array of neuron structures 
% % (each defining a different protocol) and hand them to t2n. Here we use this 
% % to simulate different current injections, but in principle this can be used 
% % for any different protocols (e.g. different channels, different synapses, different 
% % stimulation patterns etc).
% %
% amp = 0:3:15;                  % define a series of amplitudes for the electrode (here: 6 steps from 0 to 0.6 nA)
% nneuron = cell(numel(amp),1);      % initialize (clear) simulation list
% %
% %  Principally, we could copy the standard neuron structure into the neuron 
% % cell array for each protocol and then add the IClamp with different amplitude 
% % which would look like this: 
% 
% % for n = 1:numel(amp)
% %     nneuron{n} = neuron;                                                            % use neuron structure defined above as basis
% %     nneuron{n}.pp{1}.IClamp = struct('node',1,'times',[50 150],'amp',[amp(n) 0]);   % only define what is changed in for the next simulations
% % end
% %
% % However, in this case T2N would translate the complete descriptions of 
% % all mechanisms etc into hoc code for each simulation. This might be no big deal 
% % when defining only 6 simulations, but could be very time and (to lesser extent) 
% % disk space consuming when defining hundreds of  simulations. Hence, it is recommended 
% % to use the t2n_as function which tells t2n to simply reuse the hoc file of a 
% % specific simulation instance for all definitions that have not been defined 
% % in the other simulations. Here is how it is used:
% 
% nneuron{1} = neuron;                                                               % use neuron structure defined above as basis for the first simulation instance
% for n = 1:numel(amp)
%     nneuron{n}.pp{1}.IClamp = struct('node',1,'times',[50 150],'amp',[amp(n) 0]);   % now define the IClamp in all simulations (remember that all other parameters were only defined in simulation 1) and...
% end
% nneuron = t2n_as(1,nneuron);                                                        % ... use the previously defined neuron structure as basis for each other simulation
% %
% % Now T2N only writes hoc code for initialization of the IClamp point  process 
% % in each simulation and uses the hoc code about morphology,  recordings and mechanisms 
% % from the first simulation
% outs={};
% for i=1:length(tree)
%     out = t2n(nneuron,tree{i},'-w-q');                                              % execute t2n
%     outs{i}=out;
% end
% 
% %
% % Now we can plot the recorded somatic voltages for each current step
% 
% figure; hold all
% for n = 1:numel(amp)                                    % go through simulations/amplitudes
%     subplot(numel(amp),1,n)                             % make subplot for each simulation
%     hold on
%     for i=1:length(outs)
%         plot(outs{i}{n}.t,outs{i}{n}.record{1}.cell.v{1})           % plot time vs voltage of current simulation
%     end
%     hold off
%     title(sprintf('%g nA current injection',amp(n)))    % add legend
%     ylabel('Membrane pot. [mV]')
% end
% linkaxes                                                % make all axes the same size
% xlabel('Time [ms]')
% %
% % Using a "findpeaks" on the voltage vectors and plotting it against the 
% % amp vector would give us a so called f-i-relationship. However, we can save 
% % a lot of typing/formatting, if we use the the t2n functions for running current 
% % steps and plotting an f-i-relationship:
% 
% ostruct.amp = amp;%0:3:15;         % specify the same current steps [nA] as above
% ostruct.delay = 50;               % specify the same injection onset [ms] as above
% ostruct.duration = 100;           % specify the same injection duration [ms] as above
% %mkdir('tmpfolder')                % t2n_currSteps saves the result of the simulations 
%                                   % which can then be used by various plotting functions
%                                   % hence we here specify a temporary folder. this of course
%                                   % can also be a permanent folder in your simulations
% for i=1:length(tree)
%     mkdir(['tmpfolder_FIplot_',num2str(i),'/'])
%     t2n_currSteps(neuron,tree{i},['tmpfolder_FIplot_',num2str(i),'/'],ostruct); % run the current step simulations and save the results in tmpfolder
%     t2n_FIplot(['tmpfolder_FIplot_',num2str(i),'/'],neuron,ostruct); % plot an f-i-relationship from the simulation results
% end
% 
% %pause(1)
% %rmdir('tmpfolder','s')            % remove the tmpfolder
% 
% %% Tutorial E - simulation protocol: several voltage clamp steps (I-V relationship)
% % We now want to use a voltage instead of a current clamp, apply different voltage 
% % steps and plot a so called I-V curve.
% 
% amp = -140:10:-60;                 % define a series of voltage amplitude steps [mV]
% nneuron = cell(numel(amp),1);      % initialize (clear) simulation list
% %
% % We define the SEClamp electrode and again use the t2n_as function to fill 
% % up the rest of the neuron specification. Also, you see here how to record a 
% % parameter from a point process (current "i" of the SEClamp electrode in this 
% % case): Simply replace the former "cell" field (for cell recordings) with the 
% % name of the point process ("SEClamp"). The node definition specifies the node 
% % at which the point process is located in this case.
% 
% nneuron{1} = neuron;                                                               % use neuron structure defined above as basis for the first simulation instance
% nneuron{1}.record{1}.SEClamp = struct('node',1,'record','i');                      % record current "i" from SEClamp electrode which will be located at node 1
% 
% for n = 1:numel(amp)
%     % now define the SEClamp in all simulations (remember that all other parameters were only defined in simulation 1)
%     nneuron{n}.pp{1}.SEClamp = struct('node',1,'times',[0 50 150],'amp',[-80 amp(n) -80]);   
% end
% nneuron = t2n_as(1,nneuron);                                                 % then use the previously defined neuron structure as basis for each other simulation
% 
% outs={};
% for i=1:length(tree)
%     out = t2n(nneuron,tree{i},'-w-q');                                              % execute t2n
%     outs{i}=out;
% end
% 
% %
% % Now we can plot the recorded electrode currents for each voltage step
% 
% figure;
% for n = 1:numel(amp)                                    % go through simulations/amplitudes
%     subplot(ceil(sqrt(numel(amp))),ceil(sqrt(numel(amp))),n)                             % make subplot for each simulation
%     %plot(out{n}.t,out{n}.record{1}.SEClamp.i{1}*1000)        % plot time vs current of current simulation
%     hold on
%     for i=1:length(outs)
%         plot(outs{i}{n}.t,outs{i}{n}.record{1}.SEClamp.i{1}*1000)           % plot time vs voltage of current simulation
%     end
%     hold off
%     title(sprintf('%g mV VClamp',amp(n)))    % add legend
%     if rem(n-1,ceil(sqrt(numel(amp))))==0
%         ylabel('Inj. current [pA]')
%     end
%     xlim([0,200])
%     ylim([-1000 1000])
% end
% linkaxes                                                % make all axes the same size
% xlabel('Time [ms]')
% % 
% % 
% % Plotting the steady state current during the voltage step against the applied 
% % voltage step would give us the I-V relationship. However, again, we can save 
% % a lot of typing/formatting, if we use the the t2n functions for running voltage 
% % steps and plotting an I-V-relationship:
% 
% amp = -140:10:-60;         % specify the same voltage steps [mV] as above
% holding_voltage = -80;    % specify the same holding voltage [mV] as above
% duration = [50 100 50];           % specify the same injection durations [ms] as above
% %mkdir('tmpfolder')                % t2n_currSteps saves the result of the simulations 
%                                   % which can then be used by various plotting functions
%                                   % hence we here specify a temporary folder. this of course
%                                   % can also be a permanent folder in your simulations
% for i=1:length(tree)
%     mkdir(['tmpfolder_IVplot_',num2str(i),'/'])
%     t2n_voltSteps(neuron,tree{i},amp,duration,holding_voltage,['tmpfolder_IVplot_',num2str(i),'/']); % run the current step simulations and save the results in tmpfolder
%     t2n_IVplot(['tmpfolder_IVplot_',num2str(i),'/'],neuron,ostruct); % plot an f-i-relationship from the simulation results
% end
% %
% % As you can see, the I-V curve is linear at these hyperpolarized steps, 
% % as only the passive channel is conducting there.

%% Tutorial I - Synaptic stimulation, Exp2Syn synapse and a NetStim
% A more customizable synapse in NEURON is the Exp2Syn point process as it is 
% activated by another cell or point process. Hence, if we use such a synapse, 
% we have to make use of NEURON NetCon objects, which can be defined in the .con 
% field of the neuron structure. Additionally, we have to introduce some cell 
% that stimulates the synapse at given times, for which we will introduce a NetStim 
% as an artificial cell. This cell, as all real and artificial cells, has to be 
% defined in the tree structure and is quite simple to define:

vsoma = {};
time = {};
data = [];

for freq=10:10:100
    outs = {};
    synIDs = {};
    for i=1:length(tree)
        trees = {};
        trees{1} = tree{i};
        trees{2} = struct('artificial','NetStim','start',10,'interval',(1/freq)*1000,'number',ceil(freq*(neuron.params.tstop-20)/100)); % add a NetStim as an artificial cell and define the start (10 ms) the interval (15 ms) and the number (10) of spikings
        trees = t2n_writeTrees(trees,[],fullfile(pwd,['test',num2str(i),'.mtr']));                                    % tree morphologies are rewritten because this NetStim might be not written yet

        nneuron{i} = neuron;                                                         % copy standard neuron structure
        %plen = Pvec_tree(trees{1});                                                % get path length to soma at each node of morphology
        euclen = eucl_tree(trees{1});                                                % get path length to soma at each node of morphology
        %[~,synIDs{i}] = max(plen);                                                   % search for the most far away point in the morpholgy    
        %synIDs{i} = ceil(synIDs{i}/2);
        synIDs{i} = randsample(find(100<euclen<150),10);
        nneuron{i}.pp{1}.Exp2Syn = struct('node',synIDs{i},'tau1',0.2,'tau2',2.5,'e',0);% add an Exp2Syn at this location with 0.2 ms rise and 2.5 ms decay time
        nneuron{i}.record{1}.cell.node = [1,1:10:numel(tree{i}.R),synIDs{i}'];                            % record somatic voltage and voltage at synapse
        nneuron{i}.record{1}.Exp2Syn = struct('node',synIDs{i},'record','i');           % record synaptic current
        %introduce synaptic voltages
        nneuron{i}.con(1) = struct('source',struct('cell',2),'target',struct('cell',1,'pp','Exp2Syn','node',synIDs{i}),'delay',0,'threshold',0.5,'weight',.01);  % connect the NetStim (cell 2) with the target (point process Exp2Syn of cell 1 at node specified in synIDs), and add threshold/weight and delay of the connection (NetStim parameters)

        out = t2n(nneuron{i},trees,'-q');                                    % execute t2n
        outs{i} = out;
        
%         vsoma{i} = out.record{1}.cell.v{1}; 
%         time{i} = out.t;
        neuronid{i} = i;
        [~,ind] = findpeaks(out.record{1}.cell.v{1},'MinPeakHeight',0);
        freqout{i} = numel(ind)/0.3;
        freqin{i} = freq;
    end
    data = [data; [cell2mat(neuronid);cell2mat(freqin);cell2mat(freqout)]'];

    % plot the result (Vmem at soma and synapse and synaptic current)
    figure;
    subplot(2,1,1)
    hold all
    for i=1:length(outs)
        plot(outs{i}.t,outs{i}.record{1}.cell.v{1})       % plot time vs voltage at soma
        plot(outs{i}.t,outs{i}.record{1}.cell.v{synIDs{i}(1)})  % plot time vs voltage at dendrite end
    end
    legend('Soma - WT','Synapse - WT','Soma - TG','Synapse - TG')
    ylabel('Membrane potential [mV]')
    subplot(2,1,2)
    for i=1:length(outs)
        %plot(outs{i}.t,outs{i}.record{1}.Exp2Syn.i{synIDs{i}(1)})  % plot time vs synaptic current
    end
    ylabel('Synaptic current [nA]')
    xlabel('Time [ms]')

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


%% Tutorial J - Small network with Poissonian input
% Our last tutorial gives an example of how to create a small network that comprises 
% feed-forward inhibition. We use the definitions of the previous tutorial and 
% add an inhibitory integrate-and-fire point neuron that receives input from our 
% excitatory cell and and returns somatic inhibition. Furthermore, we now do not 
% use a NetStim but a VecStim, which allows us to specify a 50 Hz Poissonian spike 
% pattern that was generated with t2n_poissonSpikeGen. Last, we added a line that 
% would transform the feed-forward inhibition to a feed-back inhibition network 
% by changing the input of the inhibitory neuron from the netstim to the real 
% cell.
% 
% To speed up simulation time we furthermore enable the parallel Neuron feature 
% (cells are distributed onto different cores) + variable time step, by setting 
% nneuron.params.parallel to 4 cores and nneuron.params.cvode to 1. This feature 
% becomes more powerful if the more cells are used in the network and the more 
% cores are available to be used by T2N.

% trees = {};
% trees{1} = tree{1};
% trees{2} = struct('artificial','VecStim');                             % add a VecStim as an artificial cell. This class can have a vector played to it defining its activity
% trees{3} = struct('artificial','IntFire2','taus',15); % add a IntFire2 as an artificial cell and define its tau decay time
% trees = t2n_writeTrees(trees,[],fullfile(pwd,'test.mtr'));                                    % tree morphologies are rewritten because this NetStim might be not written yet
% 
% nneuron = neuron;                                                         % copy standard neuron structure
% nneuron.params.parallel = 4;                                              % enable parallel Neuron and set the number of used cores to 4
% nneuron.params.cvode = 1;                                                 % enable the variable time step method to fully profit from parallel Neuron (be reminded that out.t is a cell then, as time steps can be different between cells)
% plen = Pvec_tree(trees{1});                                                % get path length to soma at each node of morphology
% [~,synIDsExc] = max(plen);                                                % search for the most far away point in the morpholgy    
% nneuron.pp{1}.Exp2Syn(1) = struct('node',synIDsExc,'tau1',0.2,'tau2',2.5,'e',0);% add an Exp2Syn at this location with 0.2 ms rise and 2.5 ms decay time
% synIDsInh = 1;                                                            % choose root as point of inhibitory synapse
% nneuron.pp{1}.Exp2Syn(2) = struct('node',synIDsInh,'tau1',0.5,'tau2',5,'e',-80);% add an inhibitory (e=-80) Exp2Syn at this location with 0.2 ms rise and 2.5 ms decay time
% nneuron.record{1}.cell.node = cat(1,1,synIDsExc);                         % record somatic voltage and voltage at synapse
% nneuron.record{1}.Exp2Syn(1) = struct('node',synIDsExc,'record','i');        % record synaptic current
% nneuron.record{1}.Exp2Syn(2) = struct('node',synIDsInh,'record','i');        % record synaptic current
% nneuron.con(1) = struct('source',struct('cell',2),'target',struct('cell',1,'pp','Exp2Syn','node',synIDsExc),'delay',0.1,'threshold',0.5,'weight',0.1);  % connect the NetStim (cell 2) with the target (point process Exp2Syn of cell 1 at node specified in synIDsExc), and add threshold/weight and delay of the connection (NetStim parameters)
% nneuron.con(2) = struct('source',struct('cell',2),'target',struct('cell',3),'delay',0.1,'threshold',0.5,'weight',1);  % connect the NetStim (cell 2) with the inhibitory cell to create a feed-forward inhibitory loop, and add threshold/weight and delay of the connection (NetStim parameters)
% % nneuron.con(2) = struct('source',struct('cell',1,'node',1),'target',struct('cell',3),'delay',0,'threshold',0.5,'weight',1);  % connect the real cell with the inhibitory cell to create a feed-back inhibitory loop, and add threshold/weight and delay of the connection (NetStim parameters)
% nneuron.con(3) = struct('source',struct('cell',3),'target',struct('cell',1,'pp','Exp2Syn','node',synIDsInh),'delay',3,'threshold',0.5,'weight',0.3);  % connect the inhibitory cell with the target (point process Exp2Syn of cell 1 at node specified in synIDsInh), and add threshold/weight and delay of the connection (NetStim parameters)
% 
% [spikeMat,tvec] = t2n_poissonSpikeGen(50,neuron.params);                     % generate a 50 Hz Poissonian spike distribution
% nneuron.play{2}.cell =  struct('node',1,'play','spike','times',tvec(spikeMat));   % play this spike vector to the VecStim
% 
% out = t2n(nneuron,trees,'-w-q');                                    % execute t2n
% 
% % plot the result (Vmem at soma and synapse and synaptic current)
% fig = figure;
% subplot(4,1,2)
% hold all
% plot(out.t{1},out.record{1}.cell.v{1})       % plot time vs voltage at soma
% plot(out.t{1},out.record{1}.cell.v{synIDsExc})  % plot time vs voltage at dendrite end
% legend('Loc soma','Loc exc. synapse','Location','northoutside')
% ylabel('Membrane potential [mV]')
% subplot(4,1,3)
% plot(out.t{1},out.record{1}.Exp2Syn.i{synIDsExc})  % plot time vs synaptic current
% ylabel('Exc. syn. current [nA]')
% subplot(4,1,4)
% plot(out.t{1},out.record{1}.Exp2Syn.i{synIDsInh})  % plot time vs synaptic current
% ylabel('Inh. syn. current [nA]')
% subplot(4,1,1)
% t2n_plotRaster(spikeMat,tvec)
% title('Spikes NetStim (Poisson distr.)')
% fig.Position(4) = fig.Position(4) * 2;
% linkaxes(get(fig,'Children'),'x')