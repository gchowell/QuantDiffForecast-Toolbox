function [headline, scale, pii, rangeFrac] = quantdiffReadPII(resultFile, slots)
%QUANTDIFFREADPII  Identifiability measures of one results_W<W>_rep<r>.mat file.
%
%   [headline, scale, pii, rangeFrac] = quantdiffReadPII(resultFile, slots)
%
%   slots      slot indices into [params, X0 per fitted variable, alpha, d]
%   pii        SCI(1, slots) = log10(UB/LB) of the bootstrap 95% interval, with
%              every non-finite or complex entry set to NaN (results written
%              before November 2026 can hold Inf for a zero lower quantile and
%              complex values for a negative one)
%   rangeFrac  interval width as a fraction of the optimizer range, in [0,1]
%   scale      1 = the log-scale PII is the headline for that slot, 0 = the
%              range fraction is (lower bound <= 0); see quantdiffPII
%   headline   pii where scale == 1, rangeFrac where scale == 0
%
%   Files written by the current runners carry PII_rangeFrac and PII_scale.
%   For older files both are rebuilt from the saved interval (paramss) and
%   optimizer box (LBbox, UBbox) when those are present; otherwise rangeFrac is
%   NaN and the log scale is assumed, which reproduces the former behaviour
%   except that Inf and complex values become NaN.
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

slots = slots(:)';
info = whos('-file', resultFile);
names = {info.name};

R = load(resultFile, 'SCI');
raw = R.SCI(1, slots);
pii = real(raw);
pii(~isfinite(raw) | imag(raw) ~= 0) = NaN;

rangeFrac = nan(1, numel(slots));
scale = ones(1, numel(slots));
if all(ismember({'PII_rangeFrac', 'PII_scale'}, names))
    R2 = load(resultFile, 'PII_rangeFrac', 'PII_scale');
    rangeFrac = R2.PII_rangeFrac(1, slots);
    scale = R2.PII_scale(1, slots);
elseif all(ismember({'paramss', 'LBbox', 'UBbox'}, names))
    R2 = load(resultFile, 'paramss', 'LBbox', 'UBbox');
    c = 1 + (slots - 1) * 3;                                   % estimate column of each slot
    if max(c) + 2 <= size(R2.paramss, 2) && max(slots) <= numel(R2.LBbox)
        [~, rangeFrac, scale] = quantdiffPII(R2.paramss(1, c + 1), R2.paramss(1, c + 2), ...
            R2.LBbox(slots), R2.UBbox(slots));
    end
end
scale(~isfinite(scale)) = 1;

headline = pii;
useFrac = scale == 0;
headline(useFrac) = rangeFrac(useFrac);
end
