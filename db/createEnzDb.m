%% update in 2015
%  Enzyme definition
%  MANI, II, MGAT 1,2,3,4,5, GalT,IGNT,FucT,SiaT enzymes
%   are defined. 
residueMap                 = load('residueTypes.mat');
manResType                 = residueMap.allresidues('Man');
m3gn                       = glycanMLread('m3gn.glycoct_xml');

%% NEW DEFINITION: c1galt1  
residueMap                 = load('residueTypes.mat'); 
c1galt1                    =GTEnz([2;4;1;122]);
c1galt1.isTerminalTarget   =true;
c1galt1.resfuncgroup       =residueMap.allresidues('Gal');
c1galt1bond                =GlycanBond('1','3');
c1galt1.linkFG             = struct('anomer','b','bond',c1galt1bond); 
galnacResType              = residueMap.allresidues('GalNAc'); 
galnacBond                 =GlycanBond('?','?'); 
c1galt1.resAtt2FG          = galnacResType; 
c1galt1.linkresAtt2FG      = struct('bond', galnacBond,'anomer','a'); 
c1galt1.substNABranch      =CellArrayList;
c1galt1.substNABranch.add(glycanMLread('c1galt1subsNABranch1.glycoct_xml'));
c1galt1.substNABranch.add(glycanMLread('c1galt1subsNABranch2.glycoct_xml'));
c1galt1.substNABranch.add(glycanMLread('c1galt1subsNABranch3.glycoct_xml'));
c1galt1.substNABranch.add(glycanMLread('c1galt1subsNABranch4.glycoct_xml'));
c1galt1.substNABranch.add(glycanMLread('Ab3AN.glycoct_xml'));
c1galt1.name ='C1GALT1';


%% NEW DEFINITION: b3galt6   
residueMap                 = load('residueTypes.mat'); 
b3galt6					   =GTEnz([2;4;1;134]);
b3galt6.isTerminalTarget   =true;
b3galt6.resfuncgroup 	   =residueMap.allresidues('Gal');
b3galt6bond				   =GlycanBond('1','3');
b3galt6.linkFG             = struct('anomer','b','bond',b3galt6bond); 
galResType                 = residueMap.allresidues('Gal'); 
galnacBond                 =GlycanBond('?','?'); 
b3galt6.resAtt2FG          = galResType; 
b3galt6.linkresAtt2FG      = struct('bond', galnacBond,'anomer','b'); 
b3galt6.substNABranch      =CellArrayList;
b3galt6.substNABranch.add(glycanMLread('b3galt6substNABranch1.glycoct_xml'));
b3galt6.substNABranch.add(glycanMLread('b3galt6substNABranch2.glycoct_xml'));
b3galt6.substNABranch.add(glycanMLread('b3galt6substNABranch3.glycoct_xml'));
b3galt6.substNABranch.add(glycanMLread('b3galt6substNABranch4.glycoct_xml'));
b3galt6.substNABranch.add(glycanMLread('Ab3Ab3A.glycoct_xml'));
b3galt6.name = 'B3GALT6';

%% NEW DEFINITION: GalNAcT-A     
residueMap                 = load('residueTypes.mat'); 
galnacta                   = GTEnz([2;4;1;40]); 
galnacta.isTerminalTarget  = false; 
galnacta.resfuncgroup      = residueMap.allresidues('GalNAc'); 
galnactabond               = GlycanBond('3','1'); 
galnacta.linkFG            = struct('anomer','a','bond',galnactabond); 
galResType                 = residueMap.allresidues('Gal'); 
galBond                    = GlycanBond('?','?'); 
galnacta.resAtt2FG         =  galResType; 
galnacta.linkresAtt2FG     =  struct('bond', galBond,'anomer','b'); 
galnacta.targetNABranch    =glycanMLread('Ab4GNb2M.glycoct_xml');
galnacta.substNABranch     =glycanMLread('GalNAcT-AsubstNABranch1.glycoct_xml');
galnacta.substNABranch     =glycanMLread('GalNAcT-AsubstNABranch2.glycoct_xml');
% galnacta.substNABranch     =glycanMLread('GalNAcT-AsubstNABranch3.glycoct_xml');
galnacta.substNAResidue    =residueMap.allresidues('GalNAc');
galnacta.name ='GalNAcT-A'; 

%% NEW DEFINITION: GalT-B   
residueMap                 = load('residueTypes.mat');
galtb                      = GTEnz([2;4;1;37]);
galtb.isTerminalTarget     = false;
galtb.resfuncgroup         = residueMap.allresidues('Gal');
galtbbond                  = GlycanBond('2','1');
galtb.linkFG               = struct('anomer','a','bond',galnactabond);
galResType                 =residueMap.allresidues('Gal');
galBond                    = GlycanBond('?','?');
galtb.resAtt2FG            =  galResType;
galtb.linkresAtt2FG        =  struct('bond', galBond,'anomer','b');
galtb.targetBranch         =glycanMLread('GalTBtargetbranch.glycoct_xml');
galtb.targetNABranch       =glycanMLread('GalTBtargetNABranch1.glycoct_xml');
galtb.name = 'GalT-B'; 

%% NEW DEFINITION: FucTH1     
residueMap                 = load('residueTypes.mat');
FucTH1                     = GTEnz([2;4;1;69]);
FucTH1.isTerminalTarget    = true;
FucTH1.resfuncgroup        = residueMap.allresidues('Fuc');
FucTH1bond                 = GlycanBond('2','1');
FucTH1.linkFG              = struct('anomer','a','bond',FucTH1bond);
galResType                 = residueMap.allresidues('Gal');
galBond1				   =GlycanBond('3','1');
FucTH1.resAtt2FG           = galResType;
FucTH1.linkresAtt2FG       = struct('bond', galBond1,'anomer','b');
FucTH1.targetNABranch      =glycanMLread('Fa2Ab3GN.glycoct_xml');
FucTH1.name = 'FucTH1';

%% NEW DEFINITION: define FucTH2      
residueMap                 = load('residueTypes.mat');
FucTH2                     = GTEnz([2;4;1;69]);
FucTH2.isTerminalTarget    = true;
FucTH2.resfuncgroup        = residueMap.allresidues('Fuc');
FucTH2bond                 = GlycanBond('2','1');
FucTH2.linkFG              = struct('anomer','a','bond',FucTH2bond);
galResType                 = residueMap.allresidues('Gal');
galBond2                   =GlycanBond('4','1');
FucTH2.resAtt2FG           = galResType;
FucTH2.linkresAtt2FG       = struct('bond', galBond2,'anomer','b');
FucTH2.targetNABranch      =glycanMLread('Fa2Ab3GN.glycoct_xml');
FucTH2.name = 'FucTH2';

%% NEW DEFINITION: define FucTLe    
residueMap                 = load('residueTypes.mat');
FucTLe           		   =GTEnz([2;4;1;65]);
FucTLe.isTerminalTarget    =false;
FucTLe.resfuncgroup        =residueMap.allresidues('Fuc');
FucTLeBond                 =GlycanBond('3','1');
FucTLe.linkFG 			   =struct('anomer','a','bond',FucTLeBond);
glcnacResType			   =residueMap.allresidues('GlcNAc');
glcnacBond 				   =GlycanBond('?','?');
FucTLe.resAtt2FG 		   =glcnacResType;
FucTLe.linkresAtt2FG	   =struct('bond',glcnacBond,'anomer','b');
FucTLe.targetbranchcontain =glycanMLread('Ab4GN.glycoct_xml'); 
%FucTLe.substNAResidue     =residueMap.allresidues('Fuc');  %% maybe >1Fuc can be added    
FucTLe.substNABranch       = glycanMLread('FucTLesubstNABranch1.glycoct_xml');
FucTLe.name = 'FucTLe';

%% new define b3galt
residueMap                 = load('residueTypes.mat');
b3galt                     = GTEnz([2;4;1;62]);
b3galt.isTerminalTarget    = true;
b3galt.resfuncgroup        = residueMap.allresidues('Gal');
b3glcnacResType            = residueMap.allresidues('GlcNAc');
b3glcnacBond               = GlycanBond('?','?');
b3galt.resAtt2FG           = b3glcnacResType;
b3galt.linkresAtt2FG       = struct('bond', b3glcnacBond,'anomer','b');
b3galtbond                 = GlycanBond('3','1');
b3galt.linkFG              = struct('anomer','b','bond',b3galtbond);
b3galt.name = 'B3GalT';

%% new define IGNT
residueMap                 = load('residueTypes.mat');
manResType                 = residueMap.allresidues('Man');
IGNT                       =GTEnz([2;4;1;150]);
IGNT.isTerminalTarget      =false;
IGNT.resfuncgroup          =residueMap.allresidues('GlcNAc');
IGNTBond                   =GlycanBond('6','1');
IGNT.linkFG                =struct('anomer','b','bond',IGNTBond);
galResType                 =residueMap.allresidues('Gal');
galBond                    =GlycanBond('?','?');
IGNT.resAtt2FG             =galResType;
IGNT.linkresAtt2FG         =struct('bond',galBond,'anomer','b');
IGNT.substNABranch         =glycanMLread('IGNTsubsNABranch1.glycoct_xml');
IGNT.name = 'IGNT';

%% new define manib
residueMap                 =load('residueTypes.mat');
ManResType                 =residueMap.allresidues('Man');
manib                      =GHEnz([2;4;1;113]);
manib.isTerminalTarget     =true;
manib.resfuncgroup         =residueMap.allresidues('Man');
manib.resAtt2FG            =ManResType;
manBond                    =GlycanBond('2','1');
manunknownBond             =GlycanBond('3','1');
manib.linkresAtt2FG        =struct('anomer', 'a', 'bond', manunknownBond);
manib.linkFG               =struct('anomer','a','bond',manBond);
manib.targetBranch         =glycanMLread('Man1atargetbranch.glycoct_xml');
manib.prodMinStruct        = glycanMLread('M5.glycoct_xml'); % mannosidase I can not work on M5 glycan
manib.substMaxStruct       = glycanMLread('Man1AsubsMaxstruct.glycoct_xml');
manib.substNAResidue       =residueMap.allresidues('Glc');
manib.name = 'Man Ib';

%% new define mania
residueMap                 =load('residueTypes.mat');
ManResType                 =residueMap.allresidues('Man');
mania                      =GHEnz([2;4;1;113]);
mania.isTerminalTarget     =true;
mania.resfuncgroup         =residueMap.allresidues('Man');
mania.resAtt2FG            =ManResType;
manBond                    =GlycanBond('2','1');
manunknownBond             =GlycanBond('?','?');
mania.linkresAtt2FG        =struct('anomer', 'a', 'bond', manunknownBond);
mania.linkFG               =struct('anomer','a','bond',manBond);
mania.targetBranch         =glycanMLread('Man1atargetbranch.glycoct_xml');
mania.targetNABranch       =glycanMLread('Man1AtargetNABranch.glycoct_xml'); 
mania.prodMinStruct        = glycanMLread('M5.glycoct_xml'); % mannosidase I can not work on M5 glycan
mania.substMaxStruct       = glycanMLread('M9.glycoct_xml');
mania.substNAResidue       =residueMap.allresidues('Glc');
mania.name = 'Man I-a';

%% define sia T
siaT                       = GTEnz([2;4;99;6]);
siaT.isTerminalTarget      = true;
siaT.resfuncgroup          = residueMap.allresidues('NeuAc');
siaTbond                   = GlycanBond('3','2');
siaT.linkFG                = struct('anomer','a','bond',siaTbond);
galResType                 = residueMap.allresidues('Gal');
galBond                    = GlycanBond('4','1');
siaT.resAtt2FG             = galResType;
siaT.linkresAtt2FG         = struct('bond', galBond,'anomer','b');
siaT.substNABranch         = glycanMLread('NGlycanBisectGlcNAc.glycoct_xml');
siaT.targetbranchcontain   = glycanMLread('gntetargetbranch.glycoct_xml');
siaT.name = 'SiaT';
% enzViewer(siaT);

%% Define Fuc T
fucT8                      = GTEnz([2;4;1;68]);
fucT8.isTerminalTarget     = false;
fucT8.resfuncgroup         = residueMap.allresidues('Fuc');
fuctbond                   = GlycanBond('6','1');
fucT8.linkFG               = struct('anomer','a','bond',fuctbond);
glcnacResType              = residueMap.allresidues('GlcNAc');
glcnacBond                 = GlycanBond('?','?');
fucT8.resAtt2FG            = glcnacResType;
fucT8.linkresAtt2FG        = struct('bond', glcnacBond,'anomer','?');
fucT8.targetNABranch       = glycanMLread('NGlycanBisectGlcNAc.glycoct_xml');
fucT8.substNAResidue       = residueMap.allresidues('Gal');
fucT8.substMinStruct       = m3gn;
%  fucT8.targetBranch         = glycanMLread('fuctTargetbranch.glycoct_xml');
fucT8.substNAResidue       =residueMap.allresidues('Fuc');
fucT8.name = 'FucT8';
%enzViewer(fucT8);

%% Define MGAT1 
mgat1                      = GTEnz([2;4;1;101]);
mgat1.resfuncgroup         = residueMap.allresidues('GlcNAc');
glcnacbond                 = GlycanBond('2','1');
mgat1.linkFG               = struct('anomer','b','bond',glcnacbond);
manBond                    =  GlycanBond('3','1');
mgat1.resAtt2FG            = manResType;
mgat1.linkresAtt2FG        = struct('bond', manBond,'anomer','a');
mgat1.substMinStruct       = glycanMLread('M5.glycoct_xml');
mgat1.substMaxStruct       = glycanMLread('M5.glycoct_xml');
mgat1.targetBranch         = glycanMLread('M5_lowerbranch.glycoct_xml');
mgat1.name = 'MGAT1';
% enzViewer(mgat1);

%% Define gnte
gnte                       = GTEnz([2;4;1;149]);
gnte.isTerminalTarget      = true;
gnte.resfuncgroup          = residueMap.allresidues('GlcNAc');
galtbond                   = GlycanBond('3','1');
gnte.linkFG                = struct('anomer','b','bond',galtbond);
galResType                 = residueMap.allresidues('Gal');
galBond                    = GlycanBond('4','1');
gnte.resAtt2FG             = galResType;
gnte.linkresAtt2FG         = struct('bond', galBond,'anomer','b');
gnte.targetNABranch        = glycanMLread('NGlycanBisectGlcNAc.glycoct_xml');
gnte.name = 'GNTE';
% galt.substNABranch       =glycanMLread('NGlycanBisectGlcNAc.glycoct_xml');
%gnte.targetbranchcontain  = glycanMLread('gntetargetbranch.glycoct_xml');
% enzViewer(gnte);

%% Define MGAT2 
mgat2                      = GTEnz([2;4;1;143]);
residueMap                 = load('residueTypes.mat');
mgat2.resfuncgroup         = residueMap.allresidues('GlcNAc');
glcnacbond                 = GlycanBond('2','1');
mgat2.linkFG               = struct('anomer','b','bond',glcnacbond);
manResType                 = residueMap.allresidues('Man');
manBond                    = GlycanBond('6','1');
mgat2.resAtt2FG            = manResType;
mgat2.linkresAtt2FG        = struct('bond', manBond,'anomer','a');
m3gn                       = glycanMLread('m3gn.glycoct_xml');
mgat2.isTerminalTarget     = true;
mgat2.substMinStruct       = m3gn;
mgat2.targetBranch         = glycanMLread('mgat2actingbranch.glycoct_xml');
mgat2.substNABranch        = CellArrayList;
mgat2.substNABranch.add(glycanMLread('NGlycanBisectGlcNAc.glycoct_xml'));
mgat2.substNABranch.add(glycanMLread('mgat2substrateNAbranch.glycoct_xml'));
mgat2.name = 'MGAT2';
% enzViewer(mgat2)

%% Define MANI 
mani                       = GHEnz([3;2;1;113]);
mani.resfuncgroup          = manResType;
mani.linkFG.anomer         ='a';
manBond                    = GlycanBond('2','1');
mani.linkFG.bond           = manBond;
mani.resAtt2FG             = manResType;
manAttachBond              = GlycanBond('?','?');
mani.linkresAtt2FG         = struct('bond', manAttachBond,'anomer','?');
mani.prodMinStruct         = glycanMLread('M5.glycoct_xml'); 
mani.substMaxStruct        = glycanMLread('M9.glycoct_xml');
mani.substNAResidue        = residueMap.allresidues('Gal');
mani.name = 'Man I';
% enzViewer(mani);

%% Define MANII
manii                      = GHEnz([3;2;1;114]);
manii.resfuncgroup         = manResType;
manii.linkFG.anomer        = 'a';
manBond(1,1)               = GlycanBond('3','1');
manBond(2,1)               = GlycanBond('6','1');
manii.linkFG.bond          = manBond;
manii.resAtt2FG            = manResType;
manii.prodMinStruct        = glycanMLread('m3gn.glycoct_xml'); % mannosidase II can not work on M3 glycan
manii.substNABranch        = glycanMLread('NGlycanBisectGlcNAc.glycoct_xml');
manii.substMaxStruct       = glycanMLread('m5gn.glycoct_xml');
manii.name = 'Man II';
% enzViewer(manii)

%% Definition and visualization of MGAT3 enzyme 
mgat3                      = GTEnz([2;4;1;144]);
mgat3.resfuncgroup         = residueMap.allresidues('GlcNAc');
manResType                 = residueMap.allresidues('Man');
manBond                    = GlycanBond('4','1');
mgat3.resAtt2FG            = manResType;
mgat3.linkresAtt2FG        = struct('bond', manBond,'anomer','b');
glcnacbond                 = GlycanBond('4','1');
mgat3.linkFG               = struct('anomer','b','bond',glcnacbond);
m3gn                       = glycanMLread('m3gn.glycoct_xml');
mgat3.substMinStruct       = m3gn;
mgat3.substNABranch        = glycanMLread('NGlycanBisectGlcNAc.glycoct_xml');
mgat3.targetBranch         = glycanMLread('mgat3targetbranch.glycoct_xml');
mgat3.substNAResidue       = residueMap.allresidues('Gal');
mgat3.name = 'MGAT3';
% enzViewer(mgat3)

%% Definition and visualization of MGAT4 enzyme 
mgat4                      = GTEnz([2;4;1;145]);
mgat4.resfuncgroup         = residueMap.allresidues('GlcNAc');
manResType                 = residueMap.allresidues('Man');
manBond                    = GlycanBond('3','1');
mgat4.resAtt2FG            = manResType;
mgat4.linkresAtt2FG        = struct('bond', manBond,'anomer','a');
glcnacbond                 = GlycanBond('4','1');
mgat4.linkFG               = struct('anomer','b','bond',glcnacbond);
m3gn                       = glycanMLread('m3gn.glycoct_xml');
mgat4.substMinStruct       = m3gn;
mgat4.targetBranch         =  glycanMLread('mgat4targetbranch.glycoct_xml');
mgat4.substNABranch        = CellArrayList;
mgat4.substNABranch.add(glycanMLread('NGlycanBisectGlcNAc.glycoct_xml'));
mgat4.substNABranch.add(glycanMLread('mgat4subsNAbranch.glycoct_xml'));
mgat4.substNAResidue       = residueMap.allresidues('Gal');
mgat4.name = 'MGAT4';
%enzViewer(mgat4)

%% Definition and visualization of MGAT5 enzyme 
mgat5                      = GTEnz([2;4;1;155]);
mgat5.resfuncgroup         = residueMap.allresidues('GlcNAc');
manResType                 = residueMap.allresidues('Man');
manBond                    = GlycanBond('6','1');
mgat5.resAtt2FG            = manResType;
mgat5.linkresAtt2FG        = struct('bond', manBond,'anomer','a');
glcnacbond                 = GlycanBond('6','1');
mgat5.linkFG               = struct('anomer','b','bond',glcnacbond);
m3gn                       = glycanMLread('m3gn.glycoct_xml');
%mgat5.substMinStruct      = m3gn;
mgat5.targetBranch         = glycanMLread('mgat5targetbranch.glycoct_xml');
mgat5.substMinStruct       = glycanMLread('mgat5substMinStruct.glycoct_xml');
mgat5.substNABranch        = CellArrayList;
mgat5.substNABranch.add(glycanMLread('NGlycanBisectGlcNAc.glycoct_xml'));
mgat5.substNABranch.add(glycanMLread('mgat5substrateNAbranch.glycoct_xml'));
mgat5.substNAResidue       = residueMap.allresidues('Gal');
mgat5.name='MGAT5';

%% New define b4galt
b4galt                     = GTEnz([2;4;1;38]);
b4galt.isTerminalTarget    = true;
b4galt.resfuncgroup        = residueMap.allresidues('Gal');
glcnacResType              = residueMap.allresidues('GlcNAc');
glcnacBond                 = GlycanBond('?','1');
b4galt.resAtt2FG           = glcnacResType;
b4galt.linkresAtt2FG       = struct('bond', glcnacBond,'anomer','b');
galtbond                   = GlycanBond('4','1');
b4galt.linkFG              = struct('anomer','b','bond',galtbond);
b4galt.name='B4GalT';

ignt                       = GTEnz([2;4;1;149]);
ignt.isTerminalTarget      = true;
ignt.resfuncgroup          = residueMap.allresidues('GlcNAc');
galResType                 = residueMap.allresidues('Gal');
galBond                    = GlycanBond('4','1');
ignt.resAtt2FG             = galResType;
ignt.linkresAtt2FG         = struct(...
    'bond', galBond,'anomer','b');
glcnacbond                 = GlycanBond('3','1');
ignt.linkFG                = struct('anomer','b','bond',glcnacbond);
ignt.name='IGNT';

st3galI                    = GTEnz([2;4;99;4]);
st3galI.isTerminalTarget   = true;
st3galI.resfuncgroup       = residueMap.allresidues('NeuAc');
galResType                 = residueMap.allresidues('Gal');
galBond                    = GlycanBond('3','1');
st3galI.resAtt2FG          =  galResType;
st3galI.linkresAtt2FG      = struct( 'bond', galBond,'anomer','b');
st3bond                    = GlycanBond('3','2');
st3galI.linkFG             = struct('anomer','a','bond',st3bond);
st3galI.name = 'ST3GalI';

st3galIV                  = GTEnz([2;4;99;6]);
st3galIV.isTerminalTarget = true;
st3galIV.resfuncgroup     = residueMap.allresidues('NeuAc');
galResType                = residueMap.allresidues('Gal');
galBond                   = GlycanBond('4','1');
st3bond                   = GlycanBond('3','2');
st3galIV.linkFG           = struct('anomer','a','bond',st3bond);
st3galIV.resAtt2FG        = galResType;
st3galIV.linkresAtt2FG    = struct('bond', galBond,'anomer','b');
st3galIV.name='ST3GalIV';

% NEW st3galVI 2015
st3galVI                  = GTEnz([2;4;99;1]);
st3galVI.isTerminalTarget = true;
st3galVI.resfuncgroup     = residueMap.allresidues('NeuAc');
galResType                = residueMap.allresidues('Gal');
galBond                   = GlycanBond('4','1');
st3bond                   = GlycanBond('6','2');
st3galVI.linkFG           = struct('anomer','a','bond',st3bond);
st3galVI.resAtt2FG        = galResType;
st3galVI.linkresAtt2FG    = struct('bond', galBond,'anomer','b');
st3galVI.name = 'ST3GalVI';


fucT7                     = GTEnz([2;4;1;152]);
fucT7.isTerminalTarget    = false;
fucT7.resfuncgroup        = residueMap.allresidues('Fuc');
% fuctbond                  = GlycanBond('3','1');
fuctbond                  = GlycanBond('4','1');
fucT7.linkFG              = struct('anomer','a','bond',fuctbond);


glcnacResType             = residueMap.allresidues('GlcNAc');
glcnacBond                = GlycanBond('3','1');
fucT7.resAtt2FG           =  glcnacResType;
fucT7.linkresAtt2FG       =  struct('bond', glcnacBond,'anomer','?');
fucT7.substNAResidue      = residueMap.allresidues('Fuc');
fucT7.substMinStruct      =  glycanMLread('fucT7substmin.glycoct_xml');
fucT7.name = 'FucT7';
% enzViewer(fucT7);

enzDB= containers.Map;
enzDB('b3galt6')=b3galt6;   % new define enzyme b3galt6
enzDB('c1galt1')=c1galt1;   % new define enzyme c1galt1
enzDB('galnacta')=galnacta; % new define enzyme galnacta
enzDB('galtb')=galtb;       % new define enzyme galtb
enzDB('FucTH1')=FucTH1;     % new define enzyme FucTH1
enzDB('FucTH2')=FucTH2;     % new define enzyme FucTH2
enzDB('FucTLe')=FucTLe;     % new define enzyme FucTLe
enzDB('b3galt')=b3galt;     % new definee enzyme b3galt
enzDB('IGNT')=IGNT;         % new defined enzyme IGNT
enzDB('manib')=manib;       % new defined enzyme manib
enzDB('mania')=mania;       % new defined enzyme mania
enzDB('mgat5')=mgat5;
enzDB('mgat4')=mgat4;
enzDB('mgat3')=mgat3;
enzDB('mgat2')=mgat2;
enzDB('mgat1')=mgat1;
enzDB('manii') =manii;
enzDB('mani')  =mani;
enzDB('fucT8')=fucT8;
enzDB('siaT')=siaT;  % optional: st3galIV
enzDB('gnte')=gnte;
enzDB('b4galt')=b4galt;
enzDB('st3galI')=st3galI;
enzDB('st3galIV')=st3galIV;
enzDB('st3galVI')=st3galVI; % new 2015
enzDB('fucT7')=fucT7;

save('glyenzDB.mat','enzDB');
