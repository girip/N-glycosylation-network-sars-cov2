%% Generate array from pathway

no_reactions = length(nglycanpath.theRxns);
no_glycans = length(nglycanpath.theSpecies);

reac_array = NaN(no_reactions, 3);
no_enzs = 0;
enz_ecnos = [];
enz_names = CellArrayList;

for rxi = 1:no_reactions
    reac = nglycanpath.theRxns.get(rxi).reac.glycanStruct.name;
    prod = nglycanpath.theRxns.get(rxi).prod.glycanStruct.name;

    for ngi = 1:no_glycans
        if strcmp(reac, nglycanpath.theSpecies.get(ngi).glycanStruct.name)
            % if ~isnan(reac_array(rxi,1))
            %     disp('Multiple assignments ... error')
            % end
            reac_array(rxi,1) = ngi;
            continue
        end   
    end
    for ngi = 1:no_glycans
        if strcmp(prod, nglycanpath.theSpecies.get(ngi).glycanStruct.name)
            % if ~isnan(reac_array(rxi,2))
            %     disp('Multiple assignments ... error')
            % end
            reac_array(rxi,2) = ngi;
            continue
        end   
    end

    new_enz = 1;
    for ezi = 1:no_enzs
        if sum(nglycanpath.theRxns.get(rxi).enz.ecno' == enz_ecnos(ezi,:)) == 4 
            if ~isnan(reac_array(rxi,3))
                disp('Multiple enz assignments ... error')
            end
            %disp('Found')
            reac_array(rxi,3) = ezi;
            new_enz = 0;
        end
    end
    if new_enz || no_enzs == 0
        %disp('Adding enz')
        no_enzs = no_enzs + 1;
        reac_array(rxi,3) = no_enzs;
        enz_ecnos(no_enzs,:) = nglycanpath.theRxns.get(rxi).enz.ecno;
        enz_names.add(nglycanpath.theRxns.get(rxi).enz.name);
    end
end

%% Find the number for glycans from substrate array 

sub_gli = zeros(1,length(substrateArray))
for si = 1:length(substrateArray)
    for gli = 1:no_glycans
        if strcmp(substrateArray.get(si).name, nglycanpath.theSpecies.get(gli).name)
            sub_gli(si) = gli;
            %disp('Found !') 
            continue
        end
    end
    if sub_gli(si) == 0
        disp(strcat(['Not found ! for ' num2str(si)]))
    end
end

%%Identify initial glycans

initial_glycans = [] 
for gi = 1:no_glycans
    fgi = find(reac_array(:,2)==gi);
    if isempty(fgi)
        initial_glycans = [initial_glycans; gi];
    end
end
initial_glycans

%% Identify terminal glycans

terminal_glycans = [] 
for gi = 1:no_glycans
    fgi = find(reac_array(:,1)==gi);
    if isempty(fgi)
        terminal_glycans = [terminal_glycans; gi];
    end
end
terminal_glycans

%% Plot full network

close all 
G = digraph(reac_array(:,1),reac_array(:,2))

figure(1)
gp = plot(G,'Layout','force', 'UseGravity',true, 'linewidth',3, 'ArrowSize',8, 'NodeFontSize', 30, 'NodeFontWeight','bold')

% gp = plot(G,'Layout','circle', 'linewidth',2, 'ArrowSize',8, 'NodeFontSize', 12, 'NodeFontSize',15 , 'NodeFontWeight','bold')


highlight(gp,initial_glycans,'Marker', '*','NodeColor','k', 'MarkerSize',20,'NodeLabelColor','k')
highlight(gp,sub_gli','Marker', '.','NodeColor',[0.3, 0.3, 0.3], 'MarkerSize',20,'NodeLabelColor',[0.3, 0.3, 0.3])
% highlight(gp,terminal_glycans,'Marker', '*','NodeColor','r', 'MarkerSize',6)

for gli = 1:length(initial_glycans)
    labelnode(gp,initial_glycans(gli),num2str(initial_glycans(gli)))
end

for sgli = 1:length(sub_gli)
    if (sub_gli(sgli) ==0 )
        disp('Missing glycan pathway')
        disp(sgli)
    else
        labelnode(gp,sub_gli(sgli),num2str(sub_gli(sgli)))
    end
end

enzs = unique(reac_array(:,3));
colors = hsv(length(enzs));
% c = colormap;
% colors = c(1:fix(length(c)/length(enzs)):length(c),:);
for  enzsi = 1:length(enzs)
  rgi = find(reac_array(:,3)==enzs(enzsi));
  highlight(gp, reac_array(rgi,1),reac_array(rgi,2),'EdgeColor', colors(enzsi,:))
end

legend
hold on
for ci = 1:length(colors)
    %plot(ones(2) + ci/3,'Color',colors(ci,:))
    %text(2,1+ci/3,enz_names.get(ci),'FontSize',12)
    plot(3:4, [4+ones(1,2)] + 2*ci/3,'Color',colors(ci,:), 'linewidth',3,'linestyle','-')
    text(4,4+1+2*ci/3,enz_names.get(ci),'FontSize',12)
end
hold off

%saveas(gp,'full_network.png')
%saveas(gp,'full_network.png','png','Resolution',400)
