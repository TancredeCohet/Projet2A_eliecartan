function [paraviewPath] = getPath(paraviewCommonPaths,file_to_search,wanted_directory)
% Function retriving the path containing the file to search 
% paraviewCommonPaths is a list of cells containing usual paths where the file can be placed (in string format)
% file_to_search is the name of the file we are looking for
% wanted_directory is a label that describes well the directory 

paraviewPath = [];

% recherche du fichier dans les emplacements usuels
for i = 1:length(paraviewCommonPaths),
	if exist(strcat(paraviewCommonPaths{i},strcat('/',file_to_search)))
		paraviewPath = 	paraviewCommonPaths{i};
	end
end
% si le fichier n'est pas dans les emplacements usuels, l'utilisateur peut saisir l'emplacement du fichier si celui ci est 
% bien présent mais dans un repertoire atypique
while isempty(paraviewPath)
	question = sprintf('%s introuvable dans les emplacements usuels ; \nsaisissez nouvel emplacement : ',wanted_directory);
	answer = input(question,'s') ;
	if exist(strcat(answer,strcat('/',file_to_search)))
		paraviewPath = 	answer;
	end
end

end
