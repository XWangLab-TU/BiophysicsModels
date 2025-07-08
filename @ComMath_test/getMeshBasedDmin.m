function [D,id] = getMeshBasedDmin(coord1,coord2,Mesh,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('coord1', @(x) isnumeric(x));
ip.addRequired('coord2', @(x) isnumeric(x));
ip.addRequired('Mesh', @(x) isstruct(x));
ip.addParameter('idMesh1', [], @isnumeric);
ip.addParameter('idMesh2', [], @isnumeric);
ip.addParameter('unq1or2', 0, @isnumeric); %pick all pairs, 1: make 1 unique, 2: make 2 unique 3,: make both unique
ip.addParameter('plot_or_not', false, @islogical);
ip.parse(coord1,coord2, Mesh,varargin{:});
%--------------------------------------------------------------------------
idMesh1=ip.Results.idMesh1;
idMesh2=ip.Results.idMesh2;
unq1or2=ip.Results.unq1or2;
%--------------------------------------------------------------------------
if isempty(idMesh1)
idMesh1 = ComMath.getMeshID(coord1,Mesh.coord, Mesh.range,Mesh.d);
end
if isempty(idMesh2)
idMesh2 = ComMath.getMeshID(coord2,Mesh.coord, Mesh.range,Mesh.d);
end
%--------------------------------------------------------------------------
IdM=(idMesh1-idMesh2');
[id1,id2]=find(IdM==0);
Dpair=sqrt(sum((coord1(id1,:)-coord2(id2,:)).^2,2));
nPair=numel(Dpair);
if unq1or2==0
    id=[id1,id2];
    D=Dpair;
elseif unq1or2==1
    [Cunq,~,iCunq]=unique(id1,'row');  
    Nunq=size(Cunq,1);
    Dunq=inf(Nunq,1);    
    id=zeros(Nunq,2); 
    for iPair=1:nPair
        if Dunq(iCunq(iPair)) > Dpair(iPair)
            Dunq(iCunq(iPair))=Dpair(iPair);
            id(iCunq(iPair),:)=[id1(iPair),id2(iPair)];
        end
    end 
    D=Dunq;
elseif unq1or2==2   
    [Cunq,~,iCunq]=unique(id2,'row');  
    Nunq=size(Cunq,1);
    Dunq=inf(Nunq,1);    
    id=zeros(Nunq,2); 
    for iPair=1:nPair
        if Dunq(iCunq(iPair)) > Dpair(iPair)
            Dunq(iCunq(iPair))=Dpair(iPair);
            id(iCunq(iPair),:)=[id1(iPair),id2(iPair)];
        end
    end 
    D=Dunq;
elseif unq1or2==3
    [Cunq,~,iCunq]=unique(id1,'row');  
    Nunq=size(Cunq,1);
    Dunq=inf(Nunq,1);    
    id=zeros(Nunq,2); 
    for iPair=1:nPair
        if Dunq(iCunq(iPair)) > Dpair(iPair)
            Dunq(iCunq(iPair))=Dpair(iPair);
            id(iCunq(iPair),:)=[id1(iPair),id2(iPair)];
        end
    end 
    D=Dunq;
    id1=id(:,1);
    id2=id(:,2);
    Dpair=D;
    nPair=numel(Dpair);
    [Cunq,~,iCunq]=unique(id2,'row');  
    Nunq=size(Cunq,1);
    Dunq=inf(Nunq,1);    
    id=zeros(Nunq,2); 
    for iPair=1:nPair
        if Dunq(iCunq(iPair)) > Dpair(iPair)
            Dunq(iCunq(iPair))=Dpair(iPair);
            id(iCunq(iPair),:)=[id1(iPair),id2(iPair)];
        end
    end 
    D=Dunq;
end
%--------------------------------------------------------------------------
%%
if ip.Results.plot_or_not==true
scatter3(coord1(:,1),coord1(:,2),coord1(:,3),30,'filled','MarkerFaceColor',[1 0 0]); hold on;
scatter3(coord2(:,1),coord2(:,2),coord2(:,3),30,'filled','MarkerFaceColor',[0 0 1]); hold on;
Nunq=size(id,1);
for i=1:Nunq
plot3([coord1(id(i,1),1),coord2(id(i,2),1)],...
[coord1(id(i,1),2),coord2(id(i,2),2)],...
[coord1(id(i,1),3),coord2(id(i,2),3)],...
'linewidth',2,'color',[0 1 0]); hold on;
end
end
%--------------------------------------------------------------------------
end