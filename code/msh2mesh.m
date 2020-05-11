function msh2mesh(MshFileName)
    MshFileID = fopen(MshFileName,'r');     % Ouverture du fichier .msh (en lecture)
    if MshFileID <=0                        % Traitement du cas fichier .msh introuvable
        msg=['Le fichier de maillage : ' nomfile ' n''a pas été trouvé'];
        error(msg);
    end
    [~,FileName] = fileparts(MshFileName);      % On récupère le nom du fichier .msh pour nommer le fichier .mesh
    MeshFileName = strcat(FileName,'.mesh');    % On crée le nom complet (avec l'extension) du fichier .mesh
    MeshFileID = fopen(MeshFileName,'w');       % Ouverture (création) du fichier .mesh (en écriture)
    while ~strcmp(fgetl(MshFileID),'$Nodes'), end       % On lit successivement les lignes du fichier .msh jusqu'à arriver aux noeuds
    fprintf(MeshFileID,'MeshVersionFormatted 3\n\nDimension\n3\n\nVertices\n');     % Écriture dans le fichier .mesh les premières lignes, dimension 3 hard codée
    fprintf(MeshFileID,fgetl(MshFileID));       % Écriture dans le fichier .mesh du nombre de noeuds
    fprintf(MeshFileID,'\n');                   % Retour à la ligne pour écrire les noeuds
    current_line = fgetl(MshFileID);            % On garde la ligne courante en mémoire
    while ~strcmp(current_line,'$EndNodes')     % On parcourt tous les noeuds
        newStr = split(current_line);           % On sépare la chaîne de caractères aux espaces
        newStr = join(newStr(2:end));           % On joint les coordonnées du noeud courant (on enlève le numéro du noeud)
        newStr = string(newStr);                % cell array vers string
        fprintf(MeshFileID,newStr);             % On écrit enfin dans le fichier .mesh les coordonnées du noeud courant
        fprintf(MeshFileID,' 1\n');             % On écrit également un 1 et retour à la ligne (syntaxe .mesh)
        current_line = fgetl(MshFileID);        % On avance dans le parcourt des lignes
    end
    fprintf(MeshFileID,'\nTriangles\n');        % On écrit Triangles et retour à la ligne dans le fichier .mesh
    while ~strcmp(fgetl(MshFileID),'$Elements'), end    % On lit successivement les lignes du fichier .msh jusqu'à arriver aux éléments
    num_elements = str2double(fgetl(MshFileID));        % On garde le nombre total d'éléments en mémoire
    index_elements = 0;                         % On initialise à zéro le compteur d'éléments qui ne sont pas des triangles
    current_line = fgetl(MshFileID);            % On garde la ligne courante en mémoire
    newStr = split(current_line);               % Comme précédemment, on sépare la chaîne de caractères aux espaces
    current_element_type = newStr(2);           % On récupère le type de l'élément courant
    current_element_type = string(current_element_type);    % cell array vers string
    while ~strcmp(current_element_type,'2')     % On lit successivement les lignes du fichier .msh jusqu'à arriver aux triangles (éléments de type 2)
        current_line = fgetl(MshFileID);        % On avance dans le parcourt des lignes
        newStr = split(current_line);           % Comme précédemment, on sépare la chaîne de caractères aux espaces
        current_element_type = newStr(2);       % On récupère le type de l'élément courant
        current_element_type = string(current_element_type);    % cell array vers string
        index_elements = index_elements +1;     % On incrémente le compteur d'éléments qui ne sont pas des triangles
    end
    num_triangles = num_elements - index_elements;  % On calcule le nombre de triangles
    fprintf(MeshFileID,'%d\n',num_triangles);       % Écriture du nombre de triangles dans le fichier .mesh
    while ~strcmp(current_line,'$EndElements')      % On parcourt tous les triangles
        newStr = split(current_line);               % Comme précédemment, on sépare la chaîne de caractères aux espaces
        newStr = newStr(end-2:end);                 % On ne selectionne que les trois derniers éléments de la ligne courante (noeuds du triangle)
        newStr = join(newStr);                      % On joint les noeuds du triangle
        newStr = string(newStr);                    % cell array vers string
        fprintf(MeshFileID,newStr);                 % Écriture dans le fichier .mesh des noeuds du triangle
        fprintf(MeshFileID,' 1\n');                 % Écriture d'un 1 et retour à la ligne (syntaxe .mesh)
        current_line = fgetl(MshFileID);            % On avance dans le parcourt des lignes
    end
    fprintf(MeshFileID,'\nEnd');            % Retour à la ligne final, peut-être inutile
    fclose('all');                          % Fermture des fichiers .msh et .mesh
end
