function [F, grad_F] = free_energy(a1, b1, a2, b2, c2, lambda)
% FREE_ENERGY  Two-order-parameter toy free energy inspired by Eq. (S12).
%
% We compress the multi-q functional in Eq. (S12),
%
%   F[φ_q, X_q] = Σ_q [ (α_q + |g_q|^2/(2V_q))|X_q|^2
%                      - β_q|X_q|^4 + γ_q|X_q|^6
%                      + (g_q/V_q) φ_{-q} X_q
%                      + 1/2 (1/V_q - χ_q)|φ_q|^2 ]
%                + (quartic nonlocal term in φ) + O(|φ|^6),
%
% into a single-mode, real-scalar, local polynomial model with two fields:
%   psi1  ↔  φ   (electronic order parameter / density-like field)
%   psi2  ↔  X   (phonon/displacement-like field)
%
% Approximations / simplifications:
%  1) Keep only one representative q-mode → drop Σ_q and complex conjugation.
%  2) Ignore the φ^4 nonlocal vertex Γ_{q1q2q3q4} and higher-order φ terms.
%  3) Replace |X|^4, |X|^6 by polynomial terms in real scalar psi2.
%  4) Keep the lowest-order bilinear coupling term φ·X.
%
% Coefficient mapping (single-mode heuristic):
%   a1  ~  (1/V - χ)              (coefficient of φ^2 term)
%   a2  ~  (α + |g|^2/(2V))       (effective quadratic coefficient of X^2)
%   b2,c2 correspond to (-β, +γ)  (stabilized Landau expansion for X)
%   lambda ~ (g/V)               (bilinear coupling between φ and X)
%
% NOTE: In Eq. (S12) the X^4 term appears with a minus sign (-β|X|^4),
% while this toy code writes + (b2/4!) * psi2^4. If you want the same sign
% convention as Eq. (S12), you should pass b2 < 0 (or rewrite the term as -|b2|).
%
% Inputs:
%   a1,b1 : coefficients for psi1 (φ) polynomial (here up to quartic)
%   a2,b2,c2 : coefficients for psi2 (X) polynomial (up to 6th order for stability)
%   lambda : bilinear coupling strength between psi1 and psi2
%
% Outputs:
%   F(psi1,psi2)      : free energy function handle
%   grad_F(psi1,psi2) : gradient [∂F/∂psi1; ∂F/∂psi2]

% 定义自由能函数
F = @(psi1, psi2) (1/2)*a1*psi1.^2 + (1/factorial(4))*b1*psi1.^4 + ...
                  (1/2)*a2*psi2.^2 + (1/factorial(4))*b2*psi2.^4 + ...
                  (1/factorial(6))*c2*psi2.^6 + lambda*psi1.*psi2;

% 定义梯度函数（F 对 psi1 和 psi2 的偏导数）
grad_F = @(psi1, psi2) [a1*psi1 + (1/6)*b1*psi1.^3 + lambda*psi2; ...
                        a2*psi2 + (1/6)*b2*psi2.^3 + (1/120)*c2*psi2.^5 + lambda*psi1];
end
