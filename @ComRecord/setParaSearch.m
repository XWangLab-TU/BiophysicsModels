function [rec] = setParaSearch(rec,paraAlt,varargin)
%--------------------------------------------------------------------------
        % setParaSearch performs the preparation for @ComRecord to setup
        % parameter search, the parameter structure is stored in @ComRecord 
        % input: 
        % rec - a @ComRecord object
        % paraAlt - directory for storing the search result
        % optional:
        % see variable arguments
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%-------------------------------------------------------------------------- 
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('rec', @(x) isa(x,'ComRecord'));
ip.addRequired('paraAlt', @(x) isstruct(x));
ip.parse(rec,paraAlt,varargin{:});
%----------------------------------------------------------------------------------------
nRep=rec.nRep;
%----------------------------------------------------------------------------------------
%%
if ~isempty(paraAlt)
paraName=fieldnames(paraAlt);
nPara=numel(paraName);
nVal=zeros(nPara,1);
for i=1:nPara
    nVal(i)=numel(paraAlt.(paraName{i}));
end
nDir=1;
for i=1:nPara
    nDir=nDir*nVal(i);
end
%%
finishCount=false;
iScan=ones(1,nPara);
iTot=1;
iScanRec=iScan;
while finishCount==false
    for i=1:nPara
        if (iScan(i)<nVal(i)+1)
         iScan(i)=iScan(i)+1; 
         break;
        end
    end
        
    for i=1:nPara
        if iScan(i)==nVal(i)+1
            iScan(i)=1;
            iScan(i+1)=iScan(i+1)+1; 
%             iParaAdd=iParaAdd+1;
%             break;
        end
    end
    iScanRec=[iScanRec; iScan];
    iTot=iTot+1;

    finishCount=true;
    for i=1:nPara
        if iScan(i)~=nVal(i)
            finishCount=false;
            break;
        end
    end
end
if nDir~=iTot
    error('search range wrong!');
end
%%
dir_allNew=cell(nDir*nRep,1);

nDigit=ComMath.numDigitInt(int16(nDir));
paraDirFormat=['%0.' num2str(nDigit) 'd'];
nDigit=ComMath.numDigitInt(int16(nRep));
paraRepFormat=['%0.' num2str(nDigit) 'd'];

for i=1:nDir
   for k=1:nRep
       paraValTem=zeros(nPara,1);
       for iP=1:nPara
           paraValTem(iP)=paraAlt.(paraName{iP})(iScanRec(i,iP));
       end
   dir_allNew{(i-1)*nRep+k}=struct('dir_data',[rec.dir_all{1}.dir_data filesep 'p' num2str(i,paraDirFormat) '_rep' num2str(k,paraRepFormat)],...
                           'dir_fig',[rec.dir_all{1}.dir_fig filesep 'p' num2str(i,paraDirFormat) '_rep' num2str(k,paraRepFormat)],...
                           'dir_init',[rec.dir_all{1}.dir_init filesep 'p' num2str(i,paraDirFormat) '_rep' num2str(k,paraRepFormat)],...
                           'paraVal',paraValTem);
   end
end
else
   dir_allNew=cell(nRep,1);
   nDigit=ComMath.numDigitInt(int16(nRep));
   paraRepFormat=['%0.' num2str(nDigit) 'd'];
   for k=1:nRep
   dir_allNew{k}=struct('dir_data',[rec.dir_all{1}.dir_data filesep 'rep' num2str(k,paraRepFormat)],...
                           'dir_fig',[rec.dir_all{1}.dir_fig filesep 'rep' num2str(k,paraRepFormat)],...
                           'dir_init',[rec.dir_all{1}.dir_init filesep 'rep' num2str(k,paraRepFormat)],...
                           'paraVal',[]);
   end
end
rec.dir_all=dir_allNew;
end