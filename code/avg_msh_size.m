function [m,Nb_segments,segments,X1,X2,Y1,Y2,len_seg,coor_coll_points,coor_coll_points_pol,areas] = avg_msh_size(MeshFileName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% avg_msh_size.m:
% computes the mean element size, and other useful data in a 2D triangle
% mesh (.msh format)
% /!\ USES THE lecture_msh.m FILE /!\
%
%
% INPUT  - MeshFileName : name of the mesh file, with its .msh extension
%
% OUTPUT - m : the mean element size
%        - Nb_segments : the number of segments
%        - segments : array(nb_segments x 2), listing the segments by their
%        node relations
%        - X1,X2,Y1,Y2 : arrays of the segments coordinates
%        - len_seg : array of the segments' lengths
%        - coor_coll_points : array (number of triangles x 3) of the
%        cartesian coordinates of the collocation points (centroids)
%        - coor_coll_points_pol : array (number of triangles x 3) of the
%        polar coordinates of the collocation points (centroids)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Reading the file %s...   ',MeshFileName)
    [~,Nbtri,Coorneu,~,Numtri,~,Nbaretes,~,~]=lecture_msh(MeshFileName);
    fprintf('done.\n')
    Nb_segments = (3*Nbtri + Nbaretes)/2;   % Computes the number of segments
    segments = zeros(Nb_segments,2);        % Initializes the segments array
    index_seg = 1;                          % Keeps track of the current index in the segments array
    index_tri = 1;                          % Keeps track of the current index in the triangles array
%     if Nb_segments>10000
        
    
    f1 = waitbar(0,'1','Name','SEGMENTS (TRUE ADVANCEMENT)',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    
    setappdata(f1,'canceling',0);
    
    f2 = waitbar(0,'1','Name','TRIANGLES',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    
    setappdata(f2,'canceling',0);
    
    pos_w1=get(f1,'position');
    pos_w2=[pos_w1(1) pos_w1(2)+pos_w1(4)*1.5 pos_w1(3) pos_w1(4)];
    set(f2,'position',pos_w2,'doublebuffer','on')

    while (isequal(segments(end,:),[0 0]) && index_seg<20000)    % If the last row of segments is [0 0], it means we have not found all segments
        if getappdata(f1,'canceling')
            delete(f1)
            delete(f2)
            error('Function aborted manually.')
        end
        
        if getappdata(f2,'canceling')
            delete(f1)
            delete(f2)
            error('Function aborted manually.')
        end
        
        waitbar(index_seg/Nb_segments,f1,sprintf('%d of %d (stops at 20000) ',index_seg, Nb_segments))
        waitbar(index_tri/Nbtri,f2,sprintf('%d of %d ',index_tri, Nbtri))
        
        if ~ismember(sort(Numtri(index_tri,[1 2])), segments(1:index_seg,:), 'rows')       % Checking relation 1 2
            segments(index_seg,:) = sort(Numtri(index_tri,[1 2]));          % Working with sorted rows to avoid dupplicates
            index_seg = index_seg+1;
        end
        if ~ismember(sort(Numtri(index_tri,[1 3])), segments(1:index_seg,:), 'rows')       % Checking relation 1 3
            segments(index_seg,:) = sort(Numtri(index_tri,[1 3]));
            index_seg = index_seg+1;
        end
        if ~ismember(sort(Numtri(index_tri,[2 3])), segments(1:index_seg,:), 'rows')       % Checking relation 2 3
            segments(index_seg,:) = sort(Numtri(index_tri,[2 3]));
            index_seg = index_seg+1;
        end
        index_tri = index_tri+1;
    end
    delete(f1)
    delete(f2)
    fprintf('index_seg = %d\n', index_seg);
    
    if index_seg >= 20000
        segments = segments(1:index_seg-1,:);
        fprintf(2,'%d segments reached, stopped calculation, calculating average mesh size from the %d first segments.\n',index_seg-1,index_seg-1)
    end
    X1 = Coorneu(segments(:,1));
    X2 = Coorneu(segments(:,2));
    Y1 = Coorneu(segments(:,1),2);
    Y2 = Coorneu(segments(:,2),2);
    len_seg = sqrt((X2-X1).^2+(Y2-Y1).^2);
%     plot([X1';X2'], [Y1';Y2'])      % It is possible to plot the mesh in Matlab using these lines
%     pbaspect([1 1 1])
    m = mean(len_seg);
    
    
%     % The Coordinates and Weights of the Gauss points in the normal triangle according the numbre of the Gausse points
%     
%     coor_coll_points = zeros(Nbtri,2);% Initialization of the Cartesian Coordinates of the Gauss points in each mesh triangle
%     coor_coll_points_pol = zeros(Nbtri,2);% Initialization of the Polar Coordinates of the Gauss points in each mesh triangle
%     areas = zeros(Nbtri,1);% Initialization of the areas of each mesh triangle
%     
%     % Calculate the Cartesian Coordinates of the Gauss points in the mesh triangles
%     for t = 1:Nbtri
%         points_number = Numtri(t,:);
%         point_sommet1_coor = Coorneu(points_number(1),:);
%         point_sommet2_coor = Coorneu(points_number(2),:);
%         point_sommet3_coor = Coorneu(points_number(3),:);
%         areas(t,1) = area(point_sommet1_coor,point_sommet2_coor,point_sommet3_coor);
%         coor_coll_points(t,:) = point_sommet1_coor/3 + point_sommet2_coor/3 + point_sommet3_coor/3;
%     end
% 
%     % Transform the Cartesian Coordinates of the Gauss points in the mesh triangles to the Polar Coordinates
%     for t =1:size(coor_coll_points,1)
%         [coor_coll_points_pol(t,1),coor_coll_points_pol(t,2)] = cart2pol(coor_coll_points(t,1),coor_coll_points(t,2));
%     end
% end
coor_coll_points = 1;
coor_coll_points_pol = 2;
areas = 3;
beep
pause(1)
beep
end