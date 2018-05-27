function p = PRE_read_excel_sheet(excel_sheet,sheet,rows,numeric_fields,p_in)
%
% read excel sheet and convert into a struct
%

[num,txt,raw] = xlsread(excel_sheet,sheet);

% find field names of struct created from excel

field_names = txt(1,:);

sheet = raw(2:end,:);
nrows = size(sheet,1);
rows(rows > nrows) = [];
sheet = sheet(rows,:);

% perform eval on numeric fields (surrounded by '[' and ']')

nfields = length(field_names);
field_names_numeric = regexprep(field_names,'\d','#');
numeric_inds = find(ismember(field_names_numeric,numeric_fields));
character_inds = setdiff(1:nfields,numeric_inds);
for i = numeric_inds
    for j = 1:size(sheet,1)
        if ~isnumeric(sheet{j,i})
            sheet{j,i} = eval(['[ ' sheet{j,i} ' ]']);
        end
    end
end
for i = character_inds % turn NANs in charater fields to empty strings
    for j = 1:size(sheet,1)
        if isnan(sheet{j,i})
            sheet{j,i} = '';
        end
    end
end

% convert to struct

s = cell2struct(sheet,field_names,2);

% check which fields are numeric

for i=1:length(s)
    
    % import p_in into every element of struct
    
    if ~isempty(p_in)
        s_new(i) = import_struct(p_in,s(i));
    else
        s_new(i) = s(i);
    end
    
    % deal with session related fields
    
    session_index = 0;
    while true
        session_index = session_index+1;
        
        % find fields which end with '_<session_index>'
        
        field_inds = find(~cellfun( ...
            @isempty,regexp(field_names,['.*_' num2str(session_index)])));
        
        if(isempty(field_inds)) % if no more sessions, exit while loop
            break;
        end
        
        for field_ind = field_inds
            field_name = field_names{field_ind};
            new_field_name = regexprep(field_name,'_\d','');
            if ~isempty(s_new(i).(field_name)) && ~all(isnan(s_new(i).(field_name)))
                S(session_index).(new_field_name) = s_new(i).(field_name);
            end
        end
        
        
    end
    
    sss = s_new(i);
    if exist('S')
        sss.S = S;
        clear S;
    else
        sss.S = [];
    end
    p(i) = sss;
end

% remove fields with numeric ending

field_inds = find(~cellfun( @isempty,regexp(field_names,'.*_\d' )));
numeric_field_names = field_names(field_inds);
try 
    p = rmfield(p,numeric_field_names);
end

if ~exist('p')
    p = struct([]);
end

disp('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = import_struct(p,s_in)
%
% import struct p into struct s. Overwrite existing fields
%

s = s_in;

for field_name = fieldnames(p)';
    
    field_name = field_name{1};
    s.(field_name) = p.(field_name);
end

