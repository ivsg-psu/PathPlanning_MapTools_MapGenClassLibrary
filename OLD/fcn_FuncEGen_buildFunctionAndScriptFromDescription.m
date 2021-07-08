
%% Start building a file 
function fcn_FuncEGen_buildFunctionAndScriptFromDescription(FuncEGen)
% Check the Functions directory
if ~isfolder(fullfile(cd, 'Functions'))
    status = mkdir('Functions');
    if 0==status
        error('Unable to use the Functions directory');
    end
end

% Loop through FuncE structure, creating rebased file for each entry
for ith_function = 1:length(FuncEGen)
    
    if ith_function == 9
        disp('pause here');
    end
    
    ith_FuncE = FuncEGen(ith_function);
    
    % FUNCTION STARTS HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create file name
    fname = sprintf('fcn_%s_%s.m',ith_FuncE.class, ith_FuncE.fileNameSuffix);
    
    % Check that file does not exist
    full_fname = INTERNAL_checkIfOverwrite(fname);
    if isempty(full_fname), break, end
    
    % Open the file
    [FID_function, ~] = fopen(full_fname,'w');
        
    % Build header
    INTERNAL_fcn_writeHeader(FID_function,ith_FuncE);
    
    % Build flag area
    INTERNAL_fcn_writeFlagArea(FID_function,ith_FuncE);
    
    % Build input area
    INTERNAL_fcn_writeInputArea(FID_function,ith_FuncE);
        
    % Build main area
    INTERNAL_fcn_writeMainArea(FID_function,ith_FuncE);
    
    % Build debug area
    INTERNAL_fcn_writeDebugArea(FID_function,ith_FuncE);
    
    % Close the file
    fclose(FID_function);
    
    % SCRIPT STARTS HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create script name
    script_name = sprintf('script_test_fcn_%s_%s.m',ith_FuncE.class, ith_FuncE.fileNameSuffix);

    % Check that script does not exist
    full_script_name = INTERNAL_checkIfOverwrite(script_name);
    if isempty(full_script_name), break, end
    
    % Open the script
    [FID_script, ~] = fopen(full_script_name,'w');
    
    % Build script header
    INTERNAL_fcn_writeScriptHeader(FID_script,ith_FuncE);
    
    % Build main script area
    INTERNAL_fcn_writeScriptMainArea(FID_function,ith_FuncE);
    
    % Close the file
    fclose(FID_script);
    
end % Ends the for loop
end % Ends the function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§
%% Does file exist?
function full_fname = INTERNAL_checkIfOverwrite(fname)

full_fname = fullfile(cd, 'Functions',fname);

if isfile(full_fname) && 1==0
% if isfile(full_fname) 
    fig_warning = uifigure;
    msg = sprintf('The file: \n\n %s \n\nalready exists in directory: \n\n%s \n\nDo you want to overwrite previous versions?',fname,fullfile(cd,'Functions'));
    title = 'Confirm that you wish to overwrite.';
    selection = uiconfirm(fig_warning,...
        msg,title,...
        'Options',{'Overwrite','Cancel'},...
        'DefaultOption',2,...
        'CancelOption',2,...
        'Icon','warning');
    close(fig_warning);
    if strcmp(selection,'Cancel')
        full_fname = '';
    end
end
end

%% Write the header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _    _                _           
%  | |  | |              | |          
%  | |__| | ___  __ _  __| | ___ _ __ 
%  |  __  |/ _ \/ _` |/ _` |/ _ \ '__|
%  | |  | |  __/ (_| | (_| |  __/ |   
%  |_|  |_|\___|\__,_|\__,_|\___|_|   
%                                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                 

function INTERNAL_fcn_writeHeader(fid,FuncEdetails)
% Write start of function line
fprintf(fid,'function ');

% Write the list of outputs
if ~strcmp(FuncEdetails.Outputs(1).Name,'')
    fprintf(fid,'[ ...\n');
    
    % For each output that is required, add a variable. Create a variable
    % argument the moment one is encountered.
    Noutputs = length(FuncEdetails.Outputs);
    for ith_output = 1:Noutputs
        if FuncEdetails.Outputs(ith_output).Required == 1
            fprintf(fid,'%s',FuncEdetails.Outputs(ith_output).Name);
            if ith_output<Noutputs
                fprintf(fid,', ...\n');
            else
                fprintf(fid,' ...\n');
            end
        else
            fprintf(fid,'varargout...\n');
            break;
        end
    end
    fprintf(fid,'] = ...\n');
end

% Write the function name
fname = sprintf('fcn_%s_%s',FuncEdetails.class, FuncEdetails.fileNameSuffix);
fprintf(fid,'%s',fname);

% Write the list of inputs
if ~strcmp(FuncEdetails.Inputs(1).Name,'')
    fprintf(fid,'( ...\n');
    
    % For each input that is required, add a variable. Create a variable
    % argument the moment one is encountered.
    Ninputs = length(FuncEdetails.Inputs);
    for ith_input = 1:Ninputs
        if FuncEdetails.Inputs(ith_input).Required == 1
            fprintf(fid,'%s',FuncEdetails.Inputs(ith_input).Name);
            if ith_input<Ninputs
                fprintf(fid,', ...\n');
            else
                fprintf(fid,' ...\n');
            end
        else
            fprintf(fid,'varargin...\n');
            break;
        end
    end
    fprintf(fid,')\n');
end

% Write the commented function name
fprintf(fid,'%% %s\n',fname);

% Write the commented short description
INTERNAL_shortWrite(fid,'% ',FuncEdetails.short_description)
fprintf(fid,'%% \n');

% Write the commented long description
INTERNAL_shortWrite(fid,'% ',FuncEdetails.long_description);


% Write the format header
INTERNAL_printHeader(fid,'format');

% Write the commented list of outputs
if ~strcmp(FuncEdetails.Outputs(1).Name,'')
    fprintf(fid,'%%    [ ...\n');
    
    % For each output that is required, add a variable. Create a variable
    % argument the moment one is encountered.
    Noutputs = length(FuncEdetails.Outputs);
    for ith_output = 1:Noutputs
        if FuncEdetails.Outputs(ith_output).Required == 1
            fprintf(fid,'%%    %s',FuncEdetails.Outputs(ith_output).Name);
        else
            fprintf(fid,'%%    (%s)',FuncEdetails.Outputs(ith_output).Name);
        end
        if ith_output<Noutputs
            fprintf(fid,', ...\n');
        else
            fprintf(fid,' ...\n');
        end
    end
    fprintf(fid,'%%    ] = ...\n');
end

% Write the commented function name
fname = sprintf('fcn_%s_%s',FuncEdetails.class, FuncEdetails.fileNameSuffix);
fprintf(fid,'%%    %s',fname);

% Write the list of inputs
if ~strcmp(FuncEdetails.Inputs(1).Name,'')
    fprintf(fid,'( ...\n');
    
    % For each input that is required, add a variable. Create a variable
    % argument the moment one is encountered.
    Ninputs = length(FuncEdetails.Inputs);
    for ith_input = 1:Ninputs
        if FuncEdetails.Inputs(ith_input).Required == 1
            fprintf(fid,'%%    %s',FuncEdetails.Inputs(ith_input).Name);
        else
            fprintf(fid,'%%    (%s)',FuncEdetails.Inputs(ith_input).Name);
        end
        if ith_input<Ninputs
            fprintf(fid,', ...\n');
        else
            fprintf(fid,' ...\n');
        end
    end
    fprintf(fid,'%%    )\n');
end

% Write the inputs header
INTERNAL_printHeader(fid,'inputs');

% Write the list of input discriptions
if strcmp(FuncEdetails.Inputs(1).Name,'')
    fprintf(fid,'%%    (none)\n');
else       
    % For each input that is required, add a variable. Create a variable
    % argument the moment one is encountered.
    Ninputs = length(FuncEdetails.Inputs);
    for ith_input = 1:Ninputs
        if FuncEdetails.Inputs(ith_input).Required == 1            
            INTERNAL_shortWrite(fid,'%     ',[FuncEdetails.Inputs(ith_input).Name ': ' FuncEdetails.Inputs(ith_input).Description])                        
            fprintf(fid,'%% \n');
        else
            break
        end
    end
    
    % For each input that is NOT required, add a variable. Create a variable
    % argument the moment one is encountered.
    flag_printed_header = 0;
    for ith_input = 1:Ninputs
        if FuncEdetails.Inputs(ith_input).Required == 0
            if flag_printed_header==0
                fprintf(fid,'%%     (optional inputs)\n');
                fprintf(fid,'%%\n');
                flag_printed_header = 1;
            end
            INTERNAL_shortWrite(fid,'%     ',[FuncEdetails.Inputs(ith_input).Name ': ' FuncEdetails.Inputs(ith_input).Description])
            fprintf(fid,'%% \n');
        end
    end
end

% Write the outputs header
INTERNAL_printHeader(fid,'outputs');

% Write the list of output discriptions
if strcmp(FuncEdetails.Outputs(1).Name,'')
    fprintf(fid,'%%    (none)\n');
else       
    % For each output that is required, add a variable. Create a variable
    % argument the moment one is encountered.
    Noutputs = length(FuncEdetails.Outputs);
    for ith_output = 1:Noutputs
        if FuncEdetails.Outputs(ith_output).Required == 1            
            INTERNAL_shortWrite(fid,'%     ',[FuncEdetails.Outputs(ith_output).Name ': ' FuncEdetails.Outputs(ith_output).Description])                        
            fprintf(fid,'%% \n');
        else
            break
        end
    end
    
    % For each output that is NOT required, add a variable. Create a variable
    % argument the moment one is encountered.
    flag_printed_header = 0;
    for ith_input = 1:Noutputs
        if FuncEdetails.Outputs(ith_input).Required == 0
            if flag_printed_header==0
                fprintf(fid,'%%     (optional outputs)\n');
                fprintf(fid,'%%\n');
                flag_printed_header = 1;
            end
            INTERNAL_shortWrite(fid,'%     ',[FuncEdetails.Outputs(ith_input).Name ': ' FuncEdetails.Outputs(ith_input).Description])
            fprintf(fid,'%% \n');
        end
    end
end

% Write the dependencies header
INTERNAL_printHeader(fid,'dependencies');

% Write the list of output discriptions
if strcmp(FuncEdetails.Dependencies(1).Name,'')
    fprintf(fid,'%%    (none)\n');
else       
    % For each Dependencies, print it.
    Ndependencies = length(FuncEdetails.Dependencies);
    for ith_dependency = 1:Ndependencies
        INTERNAL_shortWrite(fid,'%     ',FuncEdetails.Dependencies(ith_dependency).Name)
        fprintf(fid,'%% \n');
    end
end


% Write the examples header
INTERNAL_printHeader(fid,'examples');

% Write the script lead
fprintf(fid,'%% See the script: script_test_%s\n',fname);
fprintf(fid,'%% for a full test suite.\n');
fprintf(fid,'%% \n');

% Write the tail
fprintf(fid,'%% This function was written on %s by %s\n',FuncEdetails.Date,FuncEdetails.Author);
fprintf(fid,'%% Questions or comments? contact %s\n',FuncEdetails.Contact);

% Break the header tail
fprintf(fid,'\n');

% Write the revision history header
INTERNAL_printHeader(fid,'revision history');
fprintf(fid,'%% %s by %s\n',FuncEdetails.Date,FuncEdetails.Author);
fprintf(fid,'%% -- first write of function\n');

% Break the header tail
fprintf(fid,'\n');

% Write the revision history header
INTERNAL_printHeader(fid,'to do');
fprintf(fid,'%% -- fill in to-do items here.\n');

% Print line break
fprintf(fid,'\n');


end

%% Write the flags area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______ _                 
%  |  ____| |                
%  | |__  | | __ _  __ _ ___ 
%  |  __| | |/ _` |/ _` / __|
%  | |    | | (_| | (_| \__ \
%  |_|    |_|\__,_|\__, |___/
%                   __/ |    
%                  |___/     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function INTERNAL_fcn_writeFlagArea(fid,FuncEdetails)
% Write start of function line
fprintf(fid,'%%%% Debugging and Input checks\n');

% List the flag option defaults, depending on flags
if FuncEdetails.Flags.CheckInputs == 1
   fprintf(fid,'flag_check_inputs = 1; %% Set equal to 1 to check the input arguments \n'); 
end
if FuncEdetails.Flags.DoPlot == 1
   fprintf(fid,'flag_do_plot = 0;      %% Set equal to 1 for plotting \n'); 
end
if FuncEdetails.Flags.DoDebug == 1
   fprintf(fid,'flag_do_debug = 0;     %% Set equal to 1 for debugging \n'); 
end

% Print line break
fprintf(fid,'\n');

% Add the debug area
if FuncEdetails.Flags.DoDebug == 1
    fprintf(fid,'if flag_do_debug\n');
    fprintf(fid,'    fig_for_debug = %.0d;\n',round(rand*1000+1));
    fprintf(fid,'    st = dbstack; %%#ok<*UNRCH>\n');
    fprintf(fid,'    fprintf(1,''STARTING function: %%s, in file: %%s\\n'',st(1).name,st(1).file);\n');
    fprintf(fid,'end \n');
end

% Print line break
fprintf(fid,'\n');

end


%% Write the inputs area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _       
%  |_   _|                 | |      
%    | |  _ __  _ __  _   _| |_ ___ 
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |                  
%              |_|                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function INTERNAL_fcn_writeInputArea(fid,FuncEdetails)

% Write the banner
writeBanner(fid,'Banner_input.txt');

% Print line break
fprintf(fid,'\n');

% Check the inputs?
if FuncEdetails.Flags.CheckInputs == 1
   fprintf(fid,'if 1 == flag_check_inputs\n');    
   
   % Print line break
   fprintf(fid,'\n');
   
   fprintf(fid,'    %% Are there the right number of inputs?\n');        
   fprintf(fid,'    if nargin < %d || nargin > %d\n',round(FuncEdetails.N_RequiredInputs),round(FuncEdetails.N_RequiredInputs+FuncEdetails.N_OptionalInputs)); 
   fprintf(fid,'        error(''Incorrect number of input arguments'')\n');
   fprintf(fid,'    end\n');

   % Print line break
   fprintf(fid,'\n');

   % For each input that is required, add a variable. Create a variable
   % argument the moment one is encountered.
   Ninputs = FuncEdetails.N_RequiredInputs;
   for ith_input = 1:Ninputs
       if isempty(FuncEdetails.Inputs(ith_input).Type)
           % do nothing
       elseif strcmp(FuncEdetails.Inputs(ith_input).Type,'char') % Char is a standard type
           fprintf(fid,'    %% Check the %s input, make sure it is characters\n',FuncEdetails.Inputs(ith_input).Name );
           fprintf(fid,'    if ~ischar(%s)\n',FuncEdetails.Inputs(ith_input).Name);
           fprintf(fid,'       error(''The variable_type_string input must be a string type, for example: ''''Path'''' '');\n');
           fprintf(fid,'    end\n');
           fprintf(fid,' \n');
       else
           fprintf(fid,'    %% Check the %s input, make sure it is ''%s'' type\n',FuncEdetails.Inputs(ith_input).Name,FuncEdetails.Inputs(ith_input).Type);
           fprintf(fid,'    fcn_%s_checkInputsToFunctions(...\n',FuncEdetails.class);
           
           % Check to see if there is a special field
           if isfield(FuncEdetails.Inputs(ith_input),'TypeOptions') && ~isempty(FuncEdetails.Inputs(ith_input).TypeOptions)
               if length(FuncEdetails.Inputs(ith_input).TypeOptions)>1
                   fprintf(fid,'        %s, ''%s'',[',FuncEdetails.Inputs(ith_input).Name,FuncEdetails.Inputs(ith_input).Type);
                   for ith_entry = 1:length(FuncEdetails.Inputs(ith_input).TypeOptions)
                       fprintf(fid,'%d ',FuncEdetails.Inputs(ith_input).TypeOptions(ith_entry));
                   end
                   fprintf(fid,']);\n');
               else
                   fprintf(fid,'        %s, ''%s'',%d);\n',FuncEdetails.Inputs(ith_input).Name,FuncEdetails.Inputs(ith_input).Type,FuncEdetails.Inputs(ith_input).TypeOptions);
               end
           else
               fprintf(fid,'        %s, ''%s'');\n',FuncEdetails.Inputs(ith_input).Name,FuncEdetails.Inputs(ith_input).Type);
           end
           fprintf(fid,' \n');
       end
   end
   
   % Close the if statement
   fprintf(fid,'end\n');
end

% Print line break
fprintf(fid,'\n');

% Check to see if user wants to show the plots
if FuncEdetails.Flags.DoPlot
    fprintf(fid,'%% Does user want to show the plots?\n');
    fprintf(fid,'if  %d== nargin\n',round(FuncEdetails.N_RequiredInputs+FuncEdetails.N_OptionalInputs));
    fprintf(fid,'    fig_num = varargin{end};\n');
    fprintf(fid,'    flag_do_plot = 1;\n');
    if FuncEdetails.Flags.DoDebug        
        fprintf(fid,'else\n');
        fprintf(fid,'    if flag_do_debug\n');
        fprintf(fid,'        fig = figure;\n');
        fprintf(fid,'        fig_for_debug = fig.Number;\n');
        fprintf(fid,'        flag_do_plot = 1;\n');
        fprintf(fid,'    end\n');
    end
    fprintf(fid,'end\n');    
end

% Print line break
fprintf(fid,'\n');

end


%% Write the main area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function INTERNAL_fcn_writeMainArea(fid,FuncEdetails) 


% Write the banner. NOTE: banner contains a special character: § which is
% used to find cut/paste area for user code in writeMain function below
writeBanner(fid,'Banner_main.txt');

% Print line break
fprintf(fid,'\n');

if ~isempty(FuncEdetails.filename_main)
    writeMain(fid,FuncEdetails.filename_main);
else
    % Tell user to put their code here
    fprintf(fid,'%% PUT YOUR CODE HERE\n');
end

% Print line break
fprintf(fid,'\n');

% Print closing character
% Print line break
fprintf(fid,'%%%s\n',char(167));

end


%% Write the debug area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _                 
%  |  __ \     | |                
%  | |  | | ___| |__  _   _  __ _ 
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function INTERNAL_fcn_writeDebugArea(fid,FuncEdetails) 

% Write the banner
writeBanner(fid,'Banner_debug.txt');

% Print line break
fprintf(fid,'\n');

% Are plots being made here?
if FuncEdetails.Flags.DoPlot  == 1
    fprintf(fid,'if flag_do_plot\n');
    fprintf(fid,'    %% Nothing to plot here\n');
    fprintf(fid,'end %% Ends the flag_do_plot if statement    \n');
end

% Print line break
fprintf(fid,'\n');

% Add the debug area?
if FuncEdetails.Flags.DoDebug == 1
    fprintf(fid,'if flag_do_debug\n');
    fprintf(fid,'    fprintf(1,''ENDING function: %%s, in file: %%s\\n\\n'',st(1).name,st(1).file);\n');
    fprintf(fid,'end\n');
end

% Print line break
fprintf(fid,'\n');

% Print line break
fprintf(fid,'\n');

% End the function!
fprintf(fid,'end %% Ends the function\n');

% Print line break
fprintf(fid,'\n');

% Write the banner for functions
writeBanner(fid,'Banner_functions.txt');

% Write the user-defined functions
writeUserFunctions(fid,FuncEdetails.filename_main);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    _____           _       _     _    _                _           
%   / ____|         (_)     | |   | |  | |              | |          
%  | (___   ___ _ __ _ _ __ | |_  | |__| | ___  __ _  __| | ___ _ __ 
%   \___ \ / __| '__| | '_ \| __| |  __  |/ _ \/ _` |/ _` |/ _ \ '__|
%   ____) | (__| |  | | |_) | |_  | |  | |  __/ (_| | (_| |  __/ |   
%  |_____/ \___|_|  |_| .__/ \__| |_|  |_|\___|\__,_|\__,_|\___|_|   
%                     | |                                            
%                     |_|                                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function INTERNAL_fcn_writeScriptHeader(fid,FuncEdetails)

% Write the header. NOTE: the last line contains a special character: §
% which is used to find cut/paste area for user code in writeScriptMain
% function below

% Write the commented script name
fprintf(fid,'%% script_test_fcn_%s_%s\n',FuncEdetails.class,FuncEdetails.fileNameSuffix);
fprintf(fid,'%% Tests: fcn_%s_%s\n',FuncEdetails.class,FuncEdetails.fileNameSuffix);


% Break the header tail
fprintf(fid,'\n');

% Write the revision history header
INTERNAL_printHeader(fid,'revision history');
fprintf(fid,'%% %s by %s\n',FuncEdetails.Date,FuncEdetails.Author);
fprintf(fid,'%% -- first write of script\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s\n',char(167));

% Break the header tail
fprintf(fid,'\n');
end

%% Write the main script area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%    _____           _       _     __  __       _       
%   / ____|         (_)     | |   |  \/  |     (_)      
%  | (___   ___ _ __ _ _ __ | |_  | \  / | __ _ _ _ __  
%   \___ \ / __| '__| | '_ \| __| | |\/| |/ _` | | '_ \ 
%   ____) | (__| |  | | |_) | |_  | |  | | (_| | | | | |
%  |_____/ \___|_|  |_| .__/ \__| |_|  |_|\__,_|_|_| |_|
%                     | |                               
%                     |_|                                                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function INTERNAL_fcn_writeScriptMainArea(fid,FuncEdetails) 

% Print line break
fprintf(fid,'\n');

if ~isempty(FuncEdetails.filename_script)
    original_name = FuncEdetails.filename_main(1:end-2); % Drop the .m at end
    new_name = sprintf('fcn_%s_%s',FuncEdetails.class,FuncEdetails.fileNameSuffix);
    writeScriptMain(fid,FuncEdetails.filename_script,original_name, new_name);
else
    % Tell user to put their code here
    fprintf(fid,'%% PUT YOUR CODE HERE\n');
end

% Print line break
fprintf(fid,'\n');

end


%%
function INTERNAL_printHeader(fid,header)
fprintf(fid,'%% \n');
fprintf(fid,'%% %s:\n',upper(header));
fprintf(fid,'%% \n');

end


%% A nice little function to constrain printing to the 75 columns total

function INTERNAL_shortWrite(fid,start_string,content)
% Writes a string that fits onto the 75 columns
% % For testing:
% start_string = '% ';
% content = '123456789012345678901234567890123456789012345678901234567890123456789012345';
content = char(content);
digits = 75 - length(start_string);

% string_format = sprintf('.{1,%.0d}',digits);
% parts = regexp(content, string_format, 'match');
% lines_to_print = parts; % Preallocate
% for ith_part = 1:length(parts)
%     lines_to_print{ith_part} = [start_string parts{ith_part}];
% end
% fprintf(fid,'%s\n',lines_to_print{:});

remainder = content;
parts = {};
while length(remainder)>digits
    current_part = remainder(1:digits);
    leftovers = remainder(digits+1:end);
    last_space = find(current_part == ' ',1,'last');
    parts = [parts, {current_part(1:last_space)}];     %#ok<AGROW>
    if last_space<length(current_part)
        remainder = [current_part(last_space+1:end) leftovers];
    else
        remainder = leftovers;
    end
end
parts = [parts, {remainder}];    

lines_to_print = parts; % Preallocate
for ith_part = 1:length(parts)
    lines_to_print{ith_part} = [start_string parts{ith_part}];
end
fprintf(fid,'%s\n',lines_to_print{:});


end



function writeBanner(fid,fname)
banner_string = fileread(fname);
fprintf(fid,'%s\n',banner_string);
end

function writeMain(fid,fname)
mainAll_string = fileread(fname);
[segments,~] = regexp(mainAll_string,char(167),'forceCellOutput','split');
fprintf(fid,'%s\n',string(segments{1}(2)));
end

function writeUserFunctions(fid,fname)
mainAll_string = fileread(fname);
[segments,~] = regexp(mainAll_string,char(167),'forceCellOutput','split');
fprintf(fid,'%s\n',string(segments{1}(4)));
end

function writeScriptMain(fid,fname,old_name,new_name)
% Grabs all characters out of the file, replaces old_name with
% new_name
mainAll_string = fileread(fname);
[segments,~] = regexp(mainAll_string,char(167),'forceCellOutput','split');
main_string = string(segments{1}(2));
main_string_fixed = regexprep(main_string,old_name,new_name);
fprintf(fid,'%s\n',main_string_fixed);
end

