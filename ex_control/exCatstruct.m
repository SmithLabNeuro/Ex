function A = exCatstruct(varargin)
% EXCATSTRUCT - concatenate structures
%
%   X = EXCATSTRUCT(S1,S2,S3,...) concates the structures S1, S2, ... into one
%   structure X.
%
%   Example:
%     A.name = 'Me' ;
%     B.income = 99999 ;
%     X = catstruct(A,B)
%     % -> X.name = 'Me' ;
%     %    X.income = 99999 ;
%
%   EXCATSTRUCT(S1,S2,'sorted') will sort the fieldnames
%   alphabetically.
%
%   EXCATSTRUCT(S1,S2,'sorted',1) will enable the error flag (so
%   that duplicate fieldnames cause an error/exit. A zero for the
%   final value, or not specifying it, will result in warnings on
%   duplicates and the last value will be used. So S2 will
%   overwrite S1, etc.
%
%   To sort the fieldnames of a structure A use:
%     A = EXCATSTRUCT(A,'sorted') ;
%
%   To concatenate two similar array of structs use simple concatenation:
%     A = dir('*.mat') ; B = dir('*.m') ; C = [A ; B] ;
%
%   When there is nothing to concatenate, the result will be an empty
%   struct (0x0 struct array with no fields).
%
%   See also CAT, STRUCT, FIELDNAMES, STRUCT2CELL

% for Matlab R13 and up
% version 2.3 (dec 2011)
% (c) Jos van der Geest
% email: jos@jasen.nl

% History
% Created:  2005
% Revisions
%   2.0 (sep 2007) removed bug when dealing with fields containing cell
%                  arrays (Thanks to Rene Willemink)
%   2.1 (sep 2008) added warning and error identifiers
%   2.2 (oct 2008) fixed error when dealing with empty structs (Thanks to
%                  Lars Barring)
%   2.3 (dec 2011) changed warning on duplicate fieldnames to error and
%                  print list of duplicates
%

N = nargin;
error(nargchk(1,Inf,N)) ;

%adding a switch option to throw an error for duplicate fields:
%-06Sep2013 ACS
if isnumeric(varargin{end})||islogical(varargin{end}), 
    errFlag = varargin{end}(1)>0;
    varargin(end) = [];
    N = N-1;
else
    errFlag = 0;
end;

if ~isstruct(varargin{end}),
    if isequal(varargin{end},'sorted'),
        sorted = 1 ;
        N = N-1 ;
        if N < 1,
            A = struct([]) ;
            return
        end
    else
        error('catstruct:InvalidArgument','Last argument should be a structure, or the string "sorted".') ;
    end
else
    sorted = 0 ;
end

FN = cell(N,1) ;
VAL = cell(N,1) ;

for ii=1:N,
    X = varargin{ii} ;
    if ~isstruct(X),
        error('catstruct:InvalidArgument',['Argument #' num2str(ii) ' is not a structure.']) ;
    end
    if ~isempty(X),
        % empty structs are ignored
        FN{ii} = fieldnames(X) ;
        VAL{ii} = struct2cell(X) ;
    end
end

FN = cat(1,FN{:}) ;
VAL = cat(1,VAL{:}) ;
[~,d] = version;
if datenum(d)<735280, %735280 is Feb 15, 2013, when R2013a was released.
    [UFN,ind] = unique(FN) ;
else
    [UFN,ind] = unique(FN,'legacy') ; %can't believe mathworks did this
end;
if numel(UFN) ~= numel(FN),    
    if ~errFlag % old method - warn on duplicates
%         warning('catstruct:DuplicatesFound','Duplicate fieldnames found. Last value is used and fields are sorted') ;
        sorted = 1 ;
    else % new method - list the duplicates and produce an error
        for ii=1:numel(UFN)
            tdup = strfind(FN,UFN{ii});
            if (length(cell2mat(tdup)) > 1)
                disp(['*** Duplicate Fieldname: ',UFN{ii}]);
            end
        end
        error('catstruct:DuplicatesFound','Duplicate fieldnames found. Exiting') ;
    end;
end

if sorted,
    VAL = VAL(ind) ;
    FN = FN(ind) ;
end

if ~isempty(FN),
    % This deals correctly with cell arrays
    A = cell2struct(VAL, FN);
else
    A = struct([]) ;
end




