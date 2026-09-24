function TGT = paper_targets()
% paper_targets -- the calibration targets, in ONE place.
%
% paper_tables.m reads them from here, so a target can never be
% updated in one report and not the other. Provenance for each is in
% Regression_moments/.
%
%   DtoY   0.6000   debt-to-output, VHLSS/SME
%   wGini  0.3720   Doan et al. (2023), 2010 value
%   exit   0.09782  Berkel-Rand-Tarp-Trifkovic (2020) Table 2.3, reproduced
%   top10  0.6321   2013 raw round, PAID labour = q101a1*(1-q101b), all firms
%   r      0.0310   Vietnam 2013 real deposit rate
%   su    -0.0269   used-capital share on log assets, 2013 SME
%   KnKu   2.5699   aggregate new-to-used capital ratio, 2013, machinery at
%                   market price. Not a fitted moment: it pins zeta through
%                   zeta = (delta+kappa)/KnKu, and the model's value is
%                   reported as a market-clearing residual.
TGT = struct('DtoY',0.6000, 'wGini',0.3720, 'exit',0.09782, ...
             'top10',0.6321, 'r',0.0310, 'su',-0.0269, 'KnKu',2.5699);
end
