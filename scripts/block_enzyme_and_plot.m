%% Blocking enzymes 

reac_array_bg = reac_array;

block_enzyme = 'FucT8'
% block_enzyme = 'Man II'
for enzi = 1:length(enz_names)
    if strcmp(enz_names.get(enzi),block_enzyme)
        bl_enzi = enzi
        continue;
    end
end

bl_reac_array = reac_array;
bl_reac_array(find(reac_array(:,3)==bl_enzi),:)=[];
bG = digraph(bl_reac_array(:,1),bl_reac_array(:,2));

bins = conncomp(bG,'Type','weak');
bin_initial = bins(initial_glycans);
nbl_reac_array = bl_reac_array(find(bins~=bin_initial),:);

% nG = rmnode(bG,find(bins~=bin_initial));

%% List glycans not included
nSubs = find(bins(sub_gli) == bin_initial);
nNotSubs = find(bins(sub_gli) ~= bin_initial);

% disp('Glycans included ..')
% for nsubi = 1:length(nSubs)
%     disp(substrateArray.get(nSubs(nsubi)).name);
% end
disp('Substrates - Formed ...')
for nsubi = 1:length(nSubs)
    disp(strcat(['Glycan :' num2str(sub_gli(nSubs(nsubi))) ' ' lf(nSubs(nsubi)).name]));
end

disp('Substrates - Not formed ...')
for nsubi = 1:length(nNotSubs)
    disp(strcat(['Glycan :' num2str(sub_gli(nNotSubs(nsubi))) ' ' lf(nNotSubs(nsubi)).name]));
    %disp(lf(nNotSubs(nsubi)).name);
end

%% Plot graph of glycan network
figure(2)
ngp = plot(bG,'Layout','force', 'UseGravity',true, 'linewidth',3, 'ArrowSize',8, 'NodeFontSize', 20, 'NodeFontWeight','bold')
title(strcat(['Blocking : ' block_enzyme]));

highlight(ngp,initial_glycans,'Marker', '*','NodeColor','k', 'MarkerSize',15,'NodeLabelColor','k')
highlight(ngp,sub_gli','Marker', '.','NodeColor',[0.3, 0.3, 0.3], 'MarkerSize',15,'NodeLabelColor',[0.3, 0.3, 0.3])
% highlight(gp,terminal_glycans,'Marker', '*','NodeColor','r', 'MarkerSize',6)

for gli = 1:length(initial_glycans)
    labelnode(ngp,initial_glycans(gli),num2str(initial_glycans(gli)))
end

for sgli = 1:length(sub_gli)
    if (sub_gli(sgli) ==0 )
        disp('Missing glycan pathway')
        disp(sgli)
    else
        labelnode(ngp,sub_gli(sgli),num2str(sub_gli(sgli)))
    end
end

enzs = unique(reac_array(:,3));
for  enzsi = 1:length(enzs)
  rgi = find(bl_reac_array(:,3)==enzs(enzsi));
  highlight(ngp, bl_reac_array(rgi,1),bl_reac_array(rgi,2),'EdgeColor', colors(enzsi,:))
end

% hold on
% for ci = 1:length(colors)
%     %plot(ones(2) + ci/3,'Color',colors(ci,:))
%     %text(2,1+ci/3,enz_names.get(ci),'FontSize',12)
%     plot(3:4, [4+ones(1,2)] + ci/3,'Color',colors(ci,:))
%     text(4,4+1+ci/3,enz_names.get(ci),'FontSize',12)
% end
% hold off 
% saveas(ngp,'blocked_network-FucT8.png')
% saveas(ngp,'blocked_network-ManII.png')
