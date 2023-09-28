function files = better_dir(folder, extension)
%BETTER_DIR Lists files in a directory and returns only files with the specified extension.
%   FILES = BETTER_DIR(FOLDER, EXTENSION) returns a structure containing the names
%   of all files in the specified directory that match the specified extension.
%   FOLDER is a string containing the path to the directory, and EXTENSION is a string
%   containing the desired file extension (e.g. 'txt', 'mat', 'csv',
%   etc.). It also works with longer file endings, say that you want to
%   return all files that end with "log_file.txt". Then simply use that in
%   place of 'extension'
%
%   The code for this function was adapted by Jonas Håkansson from a suggestion 
%   by ChatGPT, an AI language model trained by OpenAI.
%   Date 2023-02-17

% List all files in the directory
all_files = dir(folder);

% Filter for only files with the specified extension

files = all_files(endsWith({all_files.name}, extension));
end
