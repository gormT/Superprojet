function msh2mesh(MshFileName)
    MshFileID = fopen(MshFileName,'r');     % Ouverture du fichier .msh (en lecture)
    if MshFileID <=0                        % Traitement du cas fichier .msh introuvable
        msg=['Le fichier de maillage : ' nomfile ' n''a pas �t� trouv�'];
        error(msg);
    end
    [~,FileName] = fileparts(MshFileName);      % On r�cup�re le nom du fichier .msh pour nommer le fichier .mesh
    MeshFileName = strcat(FileName,'.mesh');    % On cr�e le nom complet (avec l'extension) du fichier .mesh
    MeshFileID = fopen(MeshFileName,'w');       % Ouverture (cr�ation) du fichier .mesh (en �criture)
    while ~strcmp(fgetl(MshFileID),'$Nodes'), end       % On lit successivement les lignes du fichier .msh jusqu'� arriver aux noeuds
    fprintf(MeshFileID,'MeshVersionFormatted 3\n\nDimension\n3\n\nVertices\n');     % �criture dans le fichier .mesh les premi�res lignes, dimension 3 hard cod�e
    fprintf(MeshFileID,fgetl(MshFileID));       % �criture dans le fichier .mesh du nombre de noeuds
    fprintf(MeshFileID,'\n');                   % Retour � la ligne pour �crire les noeuds
    current_line = fgetl(MshFileID);            % On garde la ligne courante en m�moire
    while ~strcmp(current_line,'$EndNodes')     % On parcourt tous les noeuds
        newStr = split(current_line);           % On s�pare la cha�ne de caract�res aux espaces
        newStr = join(newStr(2:end));           % On joint les coordonn�es du noeud courant (on enl�ve le num�ro du noeud)
        newStr = string(newStr);                % cell array vers string
        fprintf(MeshFileID,newStr);             % On �crit enfin dans le fichier .mesh les coordonn�es du noeud courant
        fprintf(MeshFileID,' 1\n');             % On �crit �galement un 1 et retour � la ligne (syntaxe .mesh)
        current_line = fgetl(MshFileID);        % On avance dans le parcourt des lignes
    end
    fprintf(MeshFileID,'\nTriangles\n');        % On �crit Triangles et retour � la ligne dans le fichier .mesh
    while ~strcmp(fgetl(MshFileID),'$Elements'), end    % On lit successivement les lignes du fichier .msh jusqu'� arriver aux �l�ments
    num_elements = str2double(fgetl(MshFileID));        % On garde le nombre total d'�l�ments en m�moire
    index_elements = 0;                         % On initialise � z�ro le compteur d'�l�ments qui ne sont pas des triangles
    current_line = fgetl(MshFileID);            % On garde la ligne courante en m�moire
    newStr = split(current_line);               % Comme pr�c�demment, on s�pare la cha�ne de caract�res aux espaces
    current_element_type = newStr(2);           % On r�cup�re le type de l'�l�ment courant
    current_element_type = string(current_element_type);    % cell array vers string
    while ~strcmp(current_element_type,'2')     % On lit successivement les lignes du fichier .msh jusqu'� arriver aux triangles (�l�ments de type 2)
        current_line = fgetl(MshFileID);        % On avance dans le parcourt des lignes
        newStr = split(current_line);           % Comme pr�c�demment, on s�pare la cha�ne de caract�res aux espaces
        current_element_type = newStr(2);       % On r�cup�re le type de l'�l�ment courant
        current_element_type = string(current_element_type);    % cell array vers string
        index_elements = index_elements +1;     % On incr�mente le compteur d'�l�ments qui ne sont pas des triangles
    end
    num_triangles = num_elements - index_elements;  % On calcule le nombre de triangles
    fprintf(MeshFileID,'%d\n',num_triangles);       % �criture du nombre de triangles dans le fichier .mesh
    while ~strcmp(current_line,'$EndElements')      % On parcourt tous les triangles
        newStr = split(current_line);               % Comme pr�c�demment, on s�pare la cha�ne de caract�res aux espaces
        newStr = newStr(end-2:end);                 % On ne selectionne que les trois derniers �l�ments de la ligne courante (noeuds du triangle)
        newStr = join(newStr);                      % On joint les noeuds du triangle
        newStr = string(newStr);                    % cell array vers string
        fprintf(MeshFileID,newStr);                 % �criture dans le fichier .mesh des noeuds du triangle
        fprintf(MeshFileID,' 1\n');                 % �criture d'un 1 et retour � la ligne (syntaxe .mesh)
        current_line = fgetl(MshFileID);            % On avance dans le parcourt des lignes
    end
    fprintf(MeshFileID,'\nEnd');            % Retour � la ligne final, peut-�tre inutile
    fclose('all');                          % Fermture des fichiers .msh et .mesh
end
