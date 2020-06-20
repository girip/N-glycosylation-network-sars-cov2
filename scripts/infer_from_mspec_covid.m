
enzdbmatfilename   = 'glyenzDB.mat';
enzdb = enzdbmatLoad(enzdbmatfilename);
mgat1 = enzdb('mgat1');
mgat2 = enzdb('mgat2');
mgat3 = enzdb('mgat3');
mgat4 = enzdb('mgat4');
mgat5 = enzdb('mgat5');

manii  = enzdb('manii');
manib    = enzdb('manib');
mania    = enzdb('mania');

gnte    = enzdb('gnte');
siaT     = enzdb('siaT');

b3galt6  = enzdb('b3galt6');
c1galt1  = enzdb('c1galt1');
galnacta = enzdb('galnacta');
galtb    = enzdb('galtb');
b3galt   = enzdb('b3galt');
IGNT     = enzdb('IGNT');
b4galt   = enzdb('b4galt');

st3galI  = enzdb('st3galI');
st3galIV = enzdb('st3galIV');
st3galVI = enzdb('st3galVI');

FucTH1   = enzdb('FucTH1');
FucTH2   = enzdb('FucTH2');
FucTLe   = enzdb('FucTLe');
fucT7    = enzdb('fucT7');
fucT8 = enzdb('fucT8');

gDispOption1    = displayset('showmass',true,'showLinkage',true,...
    'showRedEnd',true);

%  cd ../glycans/Bioarxiv_top_10_glycans/
cd ../glycans/Shajahan-glycans/

substrateArray = CellArrayList;
lf = dir('*xml')
for lfi = 1:size(lf,1)
    fname = strcat([lf(lfi).name]);
    disp(fname)
    substrateArray.add(GlycanSpecies(glycanMLread(fname)));
end

%
%% Store enzymes and glycan in CellArrayList variables
% enzArray are created to store  enzymes which might act on substrates.
enzArray = CellArrayList;
enzArray.add(manii);
enzArray.add(mgat1);
enzArray.add(mgat2);
enzArray.add(mgat3);
enzArray.add(mgat4);
enzArray.add(mgat5);
enzArray.add(fucT8);
enzArray.add(gnte);
enzArray.add(siaT);

enzArray.add(b3galt6);
enzArray.add(c1galt1);
enzArray.add(galnacta);
enzArray.add(galtb);
enzArray.add(FucTH1);
enzArray.add(FucTH2);
enzArray.add(FucTLe);
enzArray.add(b3galt);
enzArray.add(IGNT);
enzArray.add(manib);
enzArray.add(mania);
enzArray.add(b4galt);
enzArray.add(st3galI);
enzArray.add(st3galIV);
enzArray.add(st3galVI);
enzArray.add(fucT7);

enzName = CellArrayList;
enzName.add('manii');
%  enzName.add('mani');
enzName.add('mgat1');
enzName.add('mgat2');
enzName.add('mgat3');
enzName.add('mgat4');
enzName.add('mgat5');
%  enzName.add('galt');
enzName.add('fucT8');
enzName.add('gnte');
enzName.add('siaT');

enzName.add('b3galt6');
enzName.add('c1galt1');
enzName.add('galnacta');
enzName.add('galtb');
enzName.add('FucTH1');
enzName.add('FucTH2');
enzName.add('FucTLe');
enzName.add('b3galt');
enzName.add('IGNT');
enzName.add('manib');
enzName.add('mania');
enzName.add('b4galt');
enzName.add('st3galI');
enzName.add('st3galIV');
enzName.add('st3galVI');
enzName.add('fucT7');


%% Pathway reconstruction using connection inferrence
%  inferGlyConnPath command is used to construct a
%    pathway and the pathway can be visualized using glycanPathViewer
%    command.  The reconstructed network has 288 reactions and 151 species.
%    This construction takes about 5 mins to finish on Core 7 computer.
%
% [isPath,nglycanpath]=inferGlyConnPath(substrateArray, enzArray);
% glycanPathViewer(nglycanpath);

[isPath,nglycanpath]=inferGlyRevrPath(substrateArray, enzArray);
%  glycanPathViewer(nglycanpath);

