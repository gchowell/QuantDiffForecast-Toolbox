function [name, label] = quantdiffCompositeName(params)
%QUANTDIFFCOMPOSITENAME  Name of the composite (derived) quantity, safe when it is not set.
%
%   [name, label] = quantdiffCompositeName(params)
%
%   name   params.composite_name as a char row vector; '' when the options file
%          does not define the field, or leaves it empty
%   label  text for axes and legends: NAME, or 'composite parameter' when NAME
%          is empty
%
%   params.composite_name is optional. Many options files do not set it, and
%   reading the field directly stopped plotFit_ODEModel with
%   'Unrecognized field name "composite_name"'.
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

name = '';
if isstruct(params) && isfield(params, 'composite_name') && ~isempty(params.composite_name)
    v = params.composite_name;
    if isstring(v), v = char(v(1)); end
    if iscell(v) && ~isempty(v), v = v{1}; end
    if ischar(v), name = v(:)'; end
end
if isempty(name)
    label = 'composite parameter';
else
    label = name;
end
end
