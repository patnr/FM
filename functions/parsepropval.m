function prop = parsepropval(prop,varargin)
%parsepropval: Parse property/value pairs and return a structure.
%  Manages property/value pairs like MathWorks Handle Graphics functions.
%  This means that in addition to passing in Property name strings and
%  Values, you can also include structures with appropriately named fields.
%  The full, formal property names are defined by the defaults structure
%  (first input argument). This is followed by any number of property/value
%  pairs and/or structures with property names for the field names.  The
%  property name matching is case-insensitive and needs only be
%  unambiguous.
%
%  For example,
%
%    params.FileName = 'Untitled';
%    params.FileType = 'text';
%    params.Data = [1 2 3];
%    s.dat = [4 5 6];
%    parsepropval(params,'filenam','mydata.txt',s,'filety','binary')
%
%  returns a structure with the same field names as params, filled in
%  according to the property/value pairs and structures passed in.
%    ans = 
%        FileName: 'mydata.txt'
%        FileType: 'binary'
%        Data: [4 5 6]
%
%  The inputs are processed from left to right so if any property is
%  specified more than once the latest value is retained.
%
%  An error is generated if property names are ambiguous.  Values can be
%  any MATLAB variable.
%
%  Typical use is in a function with a variable number of input arguments.
%  For example,
%
%  function myfun(varargin)
%    properties.prop1 = [];
%    properties.prop2 = 'default';
%    properties = parsepropval(properties,varargin{:});


% Version: 1.0, 13 January 2009
% Author:  Douglas M. Schwarz
% Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
% Real_email = regexprep(Email,{'=','*'},{'@','.'})

% Edited by Patrick:
%		Added try-catch error handling


% Process inputs and set prop fields.
properties = fieldnames(prop);
arg_index = 1;
while arg_index <= length(varargin)
	arg = varargin{arg_index};
	if ischar(arg)
		prop_index = match_property(arg,properties);
		try
			prop.(properties{prop_index}) = varargin{arg_index + 1};
		catch exception
			error('Unmatched name/value in arguments.');
		end
		arg_index = arg_index + 2;
	elseif isstruct(arg)
		arg_fn = fieldnames(arg);
		for i = 1:length(arg_fn)
			prop_index = match_property(arg_fn{i},properties);
			prop.(properties{prop_index}) = arg.(arg_fn{i});
		end
		arg_index = arg_index + 1;
	else
		error(['Properties must be specified by property/value pairs',...
			' or structures.'])
	end
end


function prop_index = match_property(arg,properties)
prop_index = find(strcmpi(arg,properties));
if isempty(prop_index)
	prop_index = find(strncmpi(arg,properties,length(arg)));
end
if length(prop_index) ~= 1
	error('Property ''%s'' does not exist or is ambiguous.',arg)
end
