% NOTE: This function is not called by the current quasi-steady p-L pipeline.
% Reserved for a future frequency-domain p-k state-space realisation via
% rational GAF interpolation (Mayo & Antoulas 2007).  Remove this comment
% when it is wired into the stability solver.
function [L, M_loe] = loewnerInterpolation(Q_k, s_vals)
% loewnerInterpolation  Build Loewner and shifted-Loewner matrices from
%                       sampled Generalised Aerodynamic Force data.
%
%   [L, M_loe] = loewnerInterpolation(Q_k, s_vals)
%
%   Q_k    : [nModes × nModes × nK] complex GAF matrix at each sample point
%   s_vals : [1×nK] complex frequency samples  s_i = i·ω_i  (purely imaginary)
%
%   L      : Loewner matrix         [nModes·nK × nModes·nK]
%   M_loe  : Shifted-Loewner matrix [nModes·nK × nModes·nK]
%
%   The Loewner matrix interpolates Q(s) as a rational function of s,
%   enabling state-space realisation of the frequency-domain aerodynamics.
%
%   Definition (Mayo & Antoulas 2007):
%     L[i,j]     = (Q(s_i) - Q(s_j)) / (s_i - s_j)   for i ≠ j
%     M_loe[i,j] = (s_i·Q(s_i) - s_j·Q(s_j)) / (s_i - s_j)  for i ≠ j
%     Diagonal entries (i==j) are skipped (left zero).
%
%   NOTE: do NOT apply blanket small-value thresholding — it corrupts
%   legitimate small entries and breaks Hermitian structure.

nK     = length(s_vals);
nModes = size(Q_k, 1);
N      = nModes * nK;

L     = zeros(N, N);
M_loe = zeros(N, N);

for i = 1:nK
    ri = (i-1)*nModes + 1 : i*nModes;
    Qi = Q_k(:,:,i);

    for j = 1:nK
        if i == j, continue; end          % diagonal: skip (leave zero)

        ds = s_vals(i) - s_vals(j);
        if abs(ds) < 1e-14, continue; end  % guard against coincident sample points

        rj = (j-1)*nModes + 1 : j*nModes;
        Qj = Q_k(:,:,j);

        L(ri, rj)     = (Qi - Qj) / ds;
        M_loe(ri, rj) = (s_vals(i)*Qi - s_vals(j)*Qj) / ds;
    end
end
end
