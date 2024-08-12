-module(fisherl).
-export([fisher_exact/4]).

% -define(KF_GAMMA_EPS, 1.0e-14).
% -define(KF_TINY, 1.0e-290).

%-export([kf_betai_aux/3]).

%% Everyting here is ported from https://github.com/gatoravi/fisher-exact

%% Log gamma function
%% \log{\Gamma(z)}
%% AS245, 2nd algorithm, http://lib.stat.cmu.edu/apstat/245
%%
kf_lgamma(Z) ->
	X0 = 0,
	X1 = X0 + 0.1659470187408462e-06 / (Z+7),
	X2 = X1 + 0.9934937113930748e-05 / (Z+6),
	X3 = X2 - 0.1385710331296526     / (Z+5),
	X4 = X3 + 12.50734324009056      / (Z+4),
	X5 = X4 - 176.6150291498386      / (Z+3),
	X6 = X5 + 771.3234287757674      / (Z+2),
	X7 = X6 - 1259.139216722289      / (Z+1),
	X8 = X7 + 676.5203681218835      / Z,
	X9 = X8 + 0.9999999999995183,
	Result = math:log(X9) - 5.58106146679532777 - Z + (Z-0.5) * math:log(Z+6.5),
    Result.


% #define KF_GAMMA_EPS 1e-14


%     %% regularized lower incomplete gamma function, by series expansion
%     p_kf_gammap(S, Z) ->
%         Sum = 1,
%         X = 1,
%         K = 1,
%     	% for (k = 1, sum = x = 1.; k < 100; ++k) {
%     	% 	sum += (x *= z / (s + k));
%     	% 	if (x / sum < 1e-14) break;
%     	% }
%         Sum1 = p_expand(K, Sum, X, Z, S),
%         io:format("Sum1: ~f, S: ~f, Z: ~f, kflgamma: ~f~n", [Sum1, S, Z, kf_lgamma(S + 1.0)]),
%     	math:exp(S * math:log(Z) - Z - kf_lgamma(S + 1.0) + math:log(Sum1)).

%     p_expand(_K, Sum, X, _, _) when X / Sum < ?KF_GAMMA_EPS ->
%         io:format("~f~n", [Sum]),
%         Sum;
%     p_expand(K, Sum, X, Z, S) ->
%         X1 = X * Z / (S + K),
%         Sum1 = Sum + X1,
%         p_expand(K + 1, Sum1, X1, Z, S).

%     %% regularized upper incomplete gamma function, by continued fraction
%     p_kf_gammaq(S, Z) ->
%         J = 2,
%         F = 1.0 + Z - S,
%         C = F,
%         D = 0.0,
%         F1 = p_modified_lentz(J, F, C, S, Z, D, 0),
%     %	Modified Lentz's algorithm for computing continued fraction
%     %	See Numerical Recipes in C, 2nd edition, section 5.2
%     %	for (j = 1; j < 100; ++j) {
%     %		double a = j * (s - j), b = (j<<1) + 1 + z - s, d;
%     %		D = b + a * D;
%     %		if (D < KF_TINY) D = KF_TINY;
%     %		C = b + a / C;
%     %		if (C < KF_TINY) C = KF_TINY;
%     %		D = 1. / D;
%     %		d = C * D;
%     %		f *= d;
%     %		if (fabs(d - 1.) < KF_GAMMA_EPS) break;
%     %	}
%         io:format("F: ~f~n", [F1]),
%         math:exp(S * math:log(Z) - Z - kf_lgamma(S) - math:log(F1)).
%
%     p_modified_lentz(J, F, _C, _S, _Z, _D, _SmallDDiff) when J >= 100 ->
%         io:format("Quit because J limit"),
%         F;
%     p_modified_lentz(_J, F, _C, _S, _Z, _D, SmallDDiff) when SmallDDiff < ?KF_GAMMA_EPS ->
%         io:format("Quit because diff small"),
%         F;
%     p_modified_lentz(J, F, C, S, Z, D, _SmallDDiff) ->
%         A = J * (S - J),
%         B = J * 2 + 1 + Z - S,
%         D1 = B + A * D,
%         D2 = math:max(D1, ?KF_TINY),
%         C1 = B + A / C,
%         C2 = math:max(C1, ?KF_TINY),
%         D3 = 1.0 / D2,
%         SmallD = C2 * D3,
%         F1 = F * SmallD,
%         SmallDDiff = math:abs(SmallD - 1.0),
%         p_modified_lentz(J + 1, F1, C, S, Z, D3, SmallDDiff).
%
%     kf_gammap(S, Z) when (Z =< 1.0) or (Z < S) ->
%         p_kf_gammap(S, Z);
%     kf_gammap(S, Z) ->
%         1.0 - p_kf_gammaq(S, Z).
%
%
%     kf_gammaq(S, Z) when (Z =< 1.0) or (Z < S) ->
%         1.0 - p_kf_gammap(S, Z);
%     kf_gammaq(S, Z) ->
%         p_kf_gammaq(S, Z).


%     %% Regularized incomplete beta function. The method is taken from
%     %% Numerical Recipe in C, 2nd edition, section 6.4. The following web
%     %% page calculates the incomplete beta function, which equals
%     %% kf_betai(a,b,x) * gamma(a) * gamma(b) / gamma(a+b):
%     %%
%     %%   http://www.danielsoper.com/statcalc/calc36.aspx
%     %%
%     kf_betai_aux(_A, _B, X) when X == 0.0 -> 0.0;
%     kf_betai_aux(_A, _B, X) when X == 1.0 -> 1.0;
%     kf_betai_aux(A, B, X) ->
%         SmallF = 1.0,
%     	%// Modified Lentz's algorithm for computing continued fraction
%         F1 = p_lentz_continued_fraction(2, A, B, X, SmallF, SmallF, 0.0, 1.0),
%         math:exp(kf_lgamma(A + B) - kf_lgamma(A) - kf_lgamma(B) + A * math:log(X) + B * math:log(1.0 - X)) / A / F1.

%     %% Modified Lentz's algorithm for computing continued fraction
%     %%
%     %%  for (j = 1; j < 200; ++j) {
%     %%  	double aa, d;
%     %%  	int m = j>>1;
%     %%  	aa = (j&1)
%     %%  	    ? -(a + m) * (a + b + m) * x / ((a + 2*m) * (a + 2*m + 1))
%     %%  		: m * (b - m) * x / ((a + 2*m - 1) * (a + 2*m));
%     %%  	D = 1. + aa * D;
%     %%  	if (D < KF_TINY) D = KF_TINY;
%     %%  	C = 1. + aa / C;
%     %%  	if (C < KF_TINY) C = KF_TINY;
%     %%  	D = 1. / D;
%     %%  	d = C * D;
%     %%  	f *= d;
%     %%  	if (fabs(d - 1.) < KF_GAMMA_EPS) break;
%     %%  }
%     p_lentz_continued_fraction(J, _A, _B, _X, SmallF, _C, _D, _SmallDDiff) when J >= 200 ->
%         SmallF;
%     p_lentz_continued_fraction(_J, _A, _B, _X, SmallF, _C, _D, SmallDDiff) when SmallDDiff < ?KF_TINY ->
%         SmallF;
%     p_lentz_continued_fraction(J, A, B, X, SmallF, C, D, _SmallDDiff) ->
%         M = J / 2,
%         AA = p_calc_aa(J, A, B, X, M),
%         D1 = 1.0 + AA * D,
%         D2 = math:max(D1, ?KF_TINY),
%         C1 = 1.0 + AA / C,
%         C2 = math:max(C1, ?KF_TINY),
%         D3 = 1.0 / D2,
%         SmallD = C2 * D3,
%         SmallF1 = SmallF * SmallD,
%         p_lentz_continued_fraction(J + 1, A, B, X, SmallF1, C2, D3, math:abs(SmallD - 1.0)).
%
%     p_calc_aa(J, A, B, X, M) when J rem 2 == 1 ->
%     -(A + M) * (A + B + M) * X / ((A + 2*M) * (A + 2*M + 1));
%     p_calc_aa(_J, A, B, X, M) ->
%         M * (B - M) * X / ((A + 2*M - 1) * (A + 2*M)).

%% log\binom{n}{k}
lbinom(_N, K) when K == 0 ->
    0;
lbinom(N, K) when N == K ->
    0;
lbinom(N, K) ->
    kf_lgamma(N + 1) - kf_lgamma(K + 1) - kf_lgamma(N - K + 1).

%%  n11  n12  | n1_
%%  n21  n22  | n2_
%% -----------+----
%%  n_1  n_2  | n
%%
%% hypergeometric distribution
hypergeo(N11, N1_, N_1, N) ->
    math:exp(lbinom(N1_, N11) + lbinom(N - N1_, N_1 - N11) - lbinom(N, N_1)).

%% Incremental version of hypergenometric distribution.
hypergeo_acc(N11, N1_, N_1, N, {_HGN11, _HGN1_, _HGN_1, _HGN}, _HGP) when N1_ > 0;N_1 > 0;N > 0 ->
    {{N11, N1_, N_1, N}, hypergeo(N11, N1_, N_1, N)};
hypergeo_acc(N11, _N1_, _N_1, _N, {HGN11, HGN1_, HGN_1, HGN}, HGP) when ((N11 rem 11) > 0) and (N11 + HGN - HGN1_ - HGN_1 > 0) and (N11 == (HGN11 + 1)) ->
    HGP1 = HGP * (HGN1_ - HGN11) / N11
    * (HGN_1 - HGN11) / (N11 + HGN - HGN1_ - HGN_1),
    {{N11, HGN1_, HGN_1, HGN}, HGP1};
hypergeo_acc(N11, _N1_, _N_1, _N, {HGN11, HGN1_, HGN_1, HGN}, HGP) when (N11 rem 11 > 0) and (N11 + HGN - HGN1_ - HGN_1 > 0) and (N11 == (HGN11 - 1)) ->
    HGP1 = HGP * HGN11 / (HGN1_ - N11)
    * (HGN11 + HGN - HGN1_ - HGN_1) / (HGN_1 - N11),
    {{N11, HGN1_, HGN_1, HGN}, HGP1};
hypergeo_acc(N11, _N1_, _N_1, _N, {_HGN11, HGN1_, HGN_1, HGN}, _HGP) ->
    {{N11, HGN1_, HGN_1, HGN}, hypergeo(N11, HGN1_, HGN_1, HGN)}.

fisher_exact(N11, N12, N21, N22) ->
    N1_ = N11 + N12,
    N_1 = N11 + N21,
    N = N11 + N12 + N21 + N22,
    Max = min(N_1, N1_), % Max n11, for right tail % Yeah, "Max" is the lower of the two values.
    Min = N1_ + N_1 - N, %math:min(N_1, N1_), % n1_ + n_1 - n;    // not sure why n11-n22 is used instead of min(n_1,n1_)
    Min1 = max(Min, 0),% Min n11, for left tail
	if Min1 =:= Max ->
		   % No need to do test.
		   1;
	   true ->
		   % The probability of the current table.
		   {Aux, Q} = hypergeo_acc(N11, N1_, N_1, N, {0, 0, 0, 0}, 0),
		   %  Left tail.
		   {Aux1, P} = hypergeo_acc(Min1, 0, 0, 0, Aux, Q),
		   {Left, P1, I, Aux2} = sumleft(0.0, Min1 + 1, Aux1, P, Q, Max),
		   {Left2, I2} = adjustleft(P1, Q, Left, I - 1),
		   % Right tail.
		   {Aux3, PRight} = hypergeo_acc(Max, 0, 0, 0, Aux2, P1),
		   {Right, PRight1, J, _Aux4} = sumright(0.0, Max - 1, Aux3, PRight, Q),
		   {Right2, J2} = adjustright(PRight1, Q, Right, J + 1),
		   Two = gettwo(Left2, Right2),
		   {Left3, Right3} = adjust_left_right(Left2, Right2, I2, J2, N11, Q),
		   {Left3, Right3, Two}
	end.

% for (left = 0., i = min + 1; p < 0.99999999 * q && i<=max; ++i) // loop until underflow
%     left += p, p = hypergeo_acc(i, 0, 0, 0, &aux);
sumleft(Left, I, Aux, P, Q, Max) when P < 0.99999999 * Q, I =< Max ->
    {Aux1, P1} = hypergeo_acc(I, 0, 0, 0, Aux, P),
    % io:format("~.10f - i: ~p, p: ~f, aux: ~p~n", [Left, I, P, Aux]),
    sumleft(Left + P, I + 1, Aux1, P1, Q, Max);
sumleft(Left, I, Aux, P, _, _) ->
    % io:format("Final: ~.10f - i: ~p, p: ~f, aux: ~p~n", [Left, I, P, Aux]),
    {Left, P, I, Aux}.

% for (right = 0., j = max - 1; p < 0.99999999 * q && j>=0; --j) // loop until underflow
%     right += p, p = hypergeo_acc(j, 0, 0, 0, &aux);
% ++j;
sumright(Right, J, Aux, P, Q) when P < Q, J >= 0 ->
    % io:format("R: ~p, J: ~p~n", [Right, J]),
    {Aux1, P1} = hypergeo_acc(J, 0, 0, 0, Aux, P),
    sumright(Right + P, J - 1, Aux1, P1, Q);
sumright(Right, J, Aux, P, _) ->
    %io:format("Last round: R: ~p, J: ~p~n", [Right, J]),
    {Right, P, J, Aux}.

adjustleft(P, Q, Left, I) when P < 1.00000001 * Q ->
    {Left + P, I};
adjustleft(_, _, Left, I) ->
    {Left, I - 1}.


    % if (p < 1.00000001 * q) right += p;
adjustright(P, Q, Right, J) when P < 1.00000001 * Q ->
    {Right + P, J};
adjustright(_, _, Right, J) ->
    {Right, J + 1}.

gettwo(Left, Right) when Left + Right >= 1.0 ->
    1.0;
gettwo(Left, Right) ->
    Left + Right.

% // adjust left and right
% if (abs(i - n11) < abs(j - n11)) right = 1. - left + q;
% else left = 1.0 - right + q;
adjust_left_right(Left, _, I, J, N11, Q) when abs(I - N11) < abs(J - N11) ->
    Right1 = 1.0 - Left + Q,
    {Left, Right1};
adjust_left_right(_Left, Right, _, _, _, Q) ->
    {1.0 - Right + Q, Right}.


