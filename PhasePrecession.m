%% LOGAN CHIANG 2022 PHASE PRECESSION MATLAB MODEL (TSODYKS et al. 1996)

% PARAMETERS
numE = 800; numI = 200; numAll = numE + numI;

rng(0,'twister');
neuronNumber = 20;
label = linspace(0,1,numE);
%pick N neurons to observe evolution of place v phase/time
selectNeuron = randi([1 numE],1,neuronNumber); %pick which neurons?
selectLabel = label(selectNeuron); %what coordinates do those neurons have?

Vr = 0.85; thr = 1; tauM = 20; tauE = 6; tauI = 4; l = 0.08;

s_ex = 0.2; s_in = 0.7;

lmbdE = 0.03; lmbdI = 0.02; I0E = 1.02; I0I = 1.04;

J1 = 15e-3; J2 = 20e-3; sigma = 1.8;

period = 1000/8;
T = period*2; %two cycles for ease of viewing, for figs 4A/B

lapTime = 4000; dt = 1;

%SYNAPTIC WEIGHTS
syn_weights = ones(numAll,numAll)*J2;
for i = 1:numE
    for j = 1:numE
        if i<j
           syn_weights(i,j) = 1*J1*exp(-abs(label(i)-label(j))/l);
        elseif i==j
           syn_weights(i,j) = 0;
        elseif i>j
           syn_weights(i,j) = sigma*J1*exp(-abs(label(i)-label(j))/l);
        end
    end
end

tic
for travel = 1:numel(lapTime)
traversal = lapTime(travel);

external = zeros([traversal numAll]);

binary = NaN([traversal numAll]);

spike_input = zeros(numAll);

spikeCounter = zeros([traversal numAll]);

Vi = zeros([traversal numAll]); IE = zeros([traversal numAll]); II = zeros([traversal numAll]);
Vi(1,:) = Vr;

dVi = zeros([traversal-1 numAll]); dIE = zeros([traversal-1 numAll]); dII = zeros([traversal-1 numAll]);

phaseStore = NaN([traversal neuronNumber]);
labDiffStore = NaN([traversal neuronNumber]);
timeDiffStore = NaN([traversal neuronNumber]);
        
%EXTERNAL INPUT
for ext_ts = 1:traversal
    for ext = 1:numAll
        if ext <= numE
           external(ext_ts,ext) = I0E*(1+lmbdE*exp(-abs(label(ext)-ext_ts/traversal)/l));
        elseif ext>numE
           external(ext_ts,ext) = I0I*(1+lmbdI*cos(2*pi*ext_ts/period));
        end
    end
end

%BINARY RANDOM VARIABLE (PREALLOCATED FOR EFFICIENCY)
for rand_row = 1:traversal
    for rand_col = 1:numAll
        binary(rand_row,rand_col) = randi([0,1],1);
    end
end

% TIME LOOP FOR NEURONAL POTENTIAL
for timestep = 2:traversal
    
    for neuron_S = 1:numAll
        
        %SYNAPTIC CURRENT
        decayE = -IE(timestep-1,neuron_S)/tauE;
        decayI = -II(timestep-1,neuron_S)/tauI;
  
        for nloop = 1:numAll
            if nloop<=numE
            spike_input(neuron_S,nloop) = (binary(nloop,neuron_S)*s_ex)*syn_weights(nloop,neuron_S)*spikeCounter(timestep-1,nloop);
            elseif nloop>numE
            spike_input(neuron_S,nloop) = (binary(nloop,neuron_S)*s_in)*syn_weights(nloop,neuron_S)*spikeCounter(timestep-1,nloop);
            end
        end
           
        IE(timestep,neuron_S) = ((sum(spike_input(neuron_S,1:numE)))+decayE)*dt+IE(timestep-1,neuron_S);
        II(timestep,neuron_S) = ((sum(spike_input(neuron_S,numE+1:numAll)))+decayI)*dt+II(timestep-1,neuron_S);
    end
    
    Is = IE(timestep,:)-II(timestep,:);
    
    %MEMBRANE VOLTAGE
    for neuron_V = 1:numAll
        
        Vi(timestep,neuron_V) = ((-Vi(timestep-1,neuron_V)+Is(neuron_V)+external(timestep,neuron_V))/tauM)*dt+Vi(timestep-1,neuron_V);
        
        if Vi(timestep,neuron_V) >= thr
           Vi(timestep,neuron_V) = Vr;
           spikeCounter(timestep,neuron_V) = 1;
        else
           spikeCounter(timestep,neuron_V) = 0;
        end
    end
    
    %FIGURE 3
    A = spikeCounter(timestep,1:numE);
    B = any(A(:) == 1);

    if B == 0

    spiketimes = [];

    elseif B == 1

    spiketimes = find(A(:) == 1);

        xpoints1 = repmat(timestep,1,numel(spiketimes));
        ypoints1 = label(spiketimes);

        figure(1)
        x0=10;
        y0=10;
        width=810;
        height=435;
        set(gcf,'position',[x0,y0,width,height])

        hold on 
        plot(xpoints1,ypoints1,'ko','MarkerSize',3.60);
        
        ylim([0 1])
        xlim([0 1000])
        xlabel('Time (ms)');
        ylabel('Distance (m)');
        title(['Rat takes ' num2str(traversal) ' ms to traverse track at ' num2str(1000/period) ' Hz']);
        clearvars xpoints1 ypoints1

    end
    
    %FIGURE 4A/B
    for xi = 1:neuronNumber

        if spikeCounter(timestep,selectNeuron(xi)) == 1

           labDiffStore(timestep,xi) = (timestep/traversal) - selectLabel(xi);

           phaseStore(timestep,xi) = (mod(timestep,T)/T)*720;
   
        end

    end
end

hold on
figure(2)
set(gca,'FontSize',18);
stem(labDiffStore(:),phaseStore(:),'LineStyle','none','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',3);

xlim([-0.15 0.05])
xlabel('Position (m)','FontSize',20);
xticks([-0.15 -0.1 -0.05 0])
xticklabels({'-0.15','-0.1','-0.05','0'})

ylim([0 720])
ylabel('Phase','FontSize',20)
yticks([0 180 360 540 720])
yticklabels({'0','180','360','540','720'})

title([' Theta = ' num2str(1000/period) ' Hz'],'FontSize',18);

end
  
toc

%heatmap code

% bins = 36;
% nbins = [bins bins];
% [N,C] = hist3(XXX,nbins);
% contourf(C{1},C{2},N)
% 
% colormap('jet');
% 
% colorbar;

% xlim([min(XXX(:,1)) max(XXX(:,1))]);
% xlabel('Position (m)','FontSize',20);
% xticks([-0.15 -0.1 -0.05 0])
% xticklabels({'-0.15','-0.1','-0.05','0'})
% 
% ylim([min(XXX(:,2)) max(XXX(:,2))])
% ylabel('Phase','FontSize',20)
% yticks([0 180 360 540 720])
% yticklabels({'0','180','360','540','720'})
% 
% title('Theta = X Hz, Traversal Average = Y ms','FontSize',18);
