function [ver,face] = tube(obj,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('obj', @(x) isobject(x));
ip.parse(obj, varargin{:});
%--------------------------------------------------------------------------
%%
d_angle=0.04*pi;
angle=(d_angle:d_angle:2*pi)';

ver_xy=[cos(angle),sin(angle)]*10;
d=sqrt(sum((ver_xy(2,:)-ver_xy(1,:)).^2,2));
ver_xy=ver_xy/d;
d=1;
n_flat=10;
% ver_xy=[ver_xy(15:-1:5,:);[(ver_xy(5,1)+d:d:ver_xy(5,1)+10*d)' ver_xy(5,2)*ones(n_flat,1)];...
%         ver_xy(4:-1:1,:)+[ver_xy(5,1)+10*d 0];...
%         ver_xy(end:-1:16,:)+[ver_xy(5,1)+10*d 0];...
%         [(ver_xy(5,1)+10*d:-d:ver_xy(5,1)+1*d)' ver_xy(15,2)*ones(n_flat,1)]];


ver_xy_shift=[0.5*(ver_xy(1:end-1,:)+ver_xy(2:end,:));0.5*(ver_xy(1,:)+ver_xy(end,:))];
% scatter(ver_xy(:,1),ver_xy(:,2)); hold on;
% scatter(ver_xy_shift(:,1),ver_xy_shift(:,2));

ver=[];
i_turn=false;
n_ver_xy=size(ver_xy,1);
n_z=15;
for z=d:d:n_z*d
    if i_turn==false
    ver=[ver;[ver_xy ones(n_ver_xy,1)*z]];
    i_turn=true;
    else
    ver=[ver;[ver_xy_shift ones(n_ver_xy,1)*z]];
    i_turn=false;
    end
end

ver_save=ver;
ver=[ver_save(:,1),ver_save(:,3),ver_save(:,2)];

face=[];
for i_z=1:n_z-1
    i_turn=false;
    if floor(i_z/2)~=i_z/2
    for i_xy=1:n_ver_xy-1
       face=[face;[i_xy+(i_z-1)*n_ver_xy, i_xy+(i_z-1)*n_ver_xy+1, i_xy+(i_z-1)*n_ver_xy+n_ver_xy]];
    end
       face=[face;[n_ver_xy+(i_z-1)*n_ver_xy, 0+(i_z-1)*n_ver_xy+1, n_ver_xy+(i_z-1)*n_ver_xy+n_ver_xy]];
    for i_xy=1:n_ver_xy-1
       face=[face;[i_xy+(i_z-1)*n_ver_xy+1, i_xy+(i_z-1)*n_ver_xy+n_ver_xy, i_xy+(i_z-1)*n_ver_xy+n_ver_xy+1]];
    end
       face=[face;[0+(i_z-1)*n_ver_xy+1, n_ver_xy+(i_z-1)*n_ver_xy+n_ver_xy, 0+(i_z-1)*n_ver_xy+n_ver_xy+1]];
    else
    for i_xy=2:n_ver_xy
       face=[face;[i_xy+(i_z-1)*n_ver_xy, i_xy+(i_z-1)*n_ver_xy-1, i_xy+(i_z-1)*n_ver_xy+n_ver_xy]];
    end  
       face=[face;[1+(i_z-1)*n_ver_xy, n_ver_xy+(i_z-1)*n_ver_xy, 1+(i_z-1)*n_ver_xy+n_ver_xy]];
    for i_xy=2:n_ver_xy
       face=[face;[i_xy+(i_z-1)*n_ver_xy-1, i_xy+(i_z-1)*n_ver_xy+n_ver_xy-1, i_xy+(i_z-1)*n_ver_xy+n_ver_xy]];
    end 
       face=[face;[n_ver_xy+(i_z-1)*n_ver_xy, n_ver_xy+(i_z-1)*n_ver_xy+n_ver_xy, 1+(i_z-1)*n_ver_xy+n_ver_xy]];
    end
end
% patch('Vertices',ver,'Faces',face);
% scatter3(ver(:,1),ver(:,2),ver(:,3),'filled'); 

end

function [newVertices, newFaces] =  loopSubdivision(vertices, faces)
% Mesh subdivision using the Loop scheme.
%
%  Dimensions:
%    vertices: 3xnVertices
%    faces:    3xnFaces
%  
%  Author: Jesus Mena
	global edgeVertice;
    global newIndexOfVertices;
	newFaces = [];
	newVertices = vertices;
	nVertices = size(vertices,2);
	nFaces    = size(faces,2);
	edgeVertice = zeros(nVertices, nVertices, 3);
	newIndexOfVertices = nVertices;
    % ------------------------------------------------------------------------ %
	% create a matrix of edge-vertices and the new triangulation (newFaces).
    % computational complexity = O(3*nFaces)
    % 
    % * edgeVertice(x,y,1): index of the new vertice between (x,y)
    % * edgeVertice(x,y,2): index of the first opposite vertex between (x,y)
    % * edgeVertice(x,y,3): index of the second opposite vertex between (x,y)
    %
    %  0riginal vertices: va, vb, vc, vd.
    %  New vertices: vp, vq, vr.
    %
    %      vb                   vb             
    %     / \                  /  \ 
    %    /   \                vp--vq
    %   /     \              / \  / \
    % va ----- vc   ->     va-- vr --vc 
    %   \     /              \      /
    %    \   /                \    /
    %     \ /                  \  /
    %      vd                   vd               
    
	for i=1:nFaces
		[vaIndex, vbIndex, vcIndex] = deal(faces(1,i), faces(2,i), faces(3,i));
		
		vpIndex = addEdgeVertice(vaIndex, vbIndex, vcIndex);
		vqIndex = addEdgeVertice(vbIndex, vcIndex, vaIndex);
		vrIndex = addEdgeVertice(vaIndex, vcIndex, vbIndex);
		
		fourFaces = [vaIndex,vpIndex,vrIndex; vpIndex,vbIndex,vqIndex; vrIndex,vqIndex,vcIndex; vrIndex,vpIndex,vqIndex]';
		newFaces  = [newFaces, fourFaces]; 
    end;
    	
    % ------------------------------------------------------------------------ %
	% positions of the new vertices
	for v1=1:nVertices-1
		for v2=v1:nVertices
			vNIndex = edgeVertice(v1,v2,1);
            if (vNIndex~=0)
    			vNOpposite1Index = edgeVertice(v1,v2,2);
        		vNOpposite2Index = edgeVertice(v1,v2,3);
				if (vNOpposite2Index==0) % boundary case
 					newVertices(:,vNIndex) = 1/2*(vertices(:,v1)+vertices(:,v2));
				else
 					newVertices(:,vNIndex) = 3/8*(vertices(:,v1)+vertices(:,v2)) + 1/8*(vertices(:,vNOpposite1Index)+vertices(:,vNOpposite2Index));
                end;
            end;
        end;
    end;
    
	% ------------------------------------------------------------------------ %
    % adjacent vertices (using edgeVertice)
	adjVertice{nVertices} = [];
	for v=1:nVertices
		for vTmp=1:nVertices
			if (v<vTmp && edgeVertice(v,vTmp,1)~=0) || (v>vTmp && edgeVertice(vTmp,v,1)~=0)
				adjVertice{v}(end+1) = vTmp;
            end;
        end;	
    end;
    
	% ------------------------------------------------------------------------ %
    % new positions of the original vertices	
	for v=1:nVertices
		k = length(adjVertice{v});
		adjBoundaryVertices = [];
		for i=1:k
			vi = adjVertice{v}(i);
			if (vi>v) && (edgeVertice(v,vi,3)==0) || (vi<v) && (edgeVertice(vi,v,3)==0)
				adjBoundaryVertices(end+1) = vi;
			end;
		end;
		if (length(adjBoundaryVertices)==2) % boundary case
			newVertices(:,v) = 6/8*vertices(:,v) + 1/8*sum(vertices(:,adjBoundaryVertices),2);
		else
			beta = 1/k*( 5/8 - (3/8 + 1/4*cos(2*pi/k))^2 );
			newVertices(:,v) = (1-k*beta)*vertices(:,v) + beta*sum(vertices(:,(adjVertice{v})),2); 
		end;
    end;
 	
end
% ---------------------------------------------------------------------------- %
function vNIndex = addEdgeVertice(v1Index, v2Index, v3Index)
	global edgeVertice;
	global newIndexOfVertices;
	if (v1Index>v2Index) % setting: v1 <= v2
		vTmp = v1Index;
		v1Index = v2Index;
		v2Index = vTmp;
	end;
	
	if (edgeVertice(v1Index, v2Index, 1)==0)  % new vertex
		newIndexOfVertices = newIndexOfVertices+1;
		edgeVertice(v1Index, v2Index, 1) = newIndexOfVertices;
		edgeVertice(v1Index, v2Index, 2) = v3Index;
	else
		edgeVertice(v1Index, v2Index, 3) = v3Index;
	end;
	vNIndex = edgeVertice(v1Index, v2Index, 1);
    return;
end