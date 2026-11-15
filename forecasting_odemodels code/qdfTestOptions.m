function qdfTestOptions(src, dst, varargin)
%QDFTESTOPTIONS  Make a modified copy of an options file for a quick test run.
%
%   qdfTestOptions(src, dst, old1, new1, old2, new2, ...)
%
%   Copies the options file SRC (name without .m) to DST.m, renames the
%   function inside it to DST, and replaces each text OLDk by NEWk. It stops
%   with an error if an OLDk text is not found, so a test never runs with an
%   unchanged setting by mistake.
%
%   Example: a fast copy of the negative binomial example (20 bootstrap fits)
%     qdfTestOptions('options_fit_SEIR_flu1918_dist1_3', 'options_test_fit', ...
%                    'B = 300;', 'B = 20;');
%
%   Delete the DST.m files when the tests are done.

if mod(numel(varargin), 2) ~= 0
    error('qdfTestOptions:pairs', 'Replacements must come in pairs: old text, new text.');
end

txt = fileread([src '.m']);
if isempty(strfind(txt, src)) %#ok<STREMP>
    error('qdfTestOptions:name', 'The function name %s was not found inside %s.m.', src, src);
end
txt = strrep(txt, src, dst);

for k = 1:2:numel(varargin)
    oldText = varargin{k};
    newText = varargin{k+1};
    if isempty(strfind(txt, oldText)) %#ok<STREMP>
        error('qdfTestOptions:notFound', 'Text not found in %s.m: "%s"', src, oldText);
    end
    txt = strrep(txt, oldText, newText);
end

fid = fopen([dst '.m'], 'w');
if fid < 0
    error('qdfTestOptions:write', 'Cannot write %s.m in the current folder.', dst);
end
fwrite(fid, txt);
fclose(fid);
rehash
fprintf('Created %s.m from %s.m\n', dst, src);
end
