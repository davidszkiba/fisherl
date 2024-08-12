-module(fisherl).
-export([fisher_exact/4]).

%% Everyting here is ported from https://github.com/gatoravi/fisher-exact

%% Log gamma function
%%
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

%% log\binom{n}{k}
lbinom(_N, K) when K == 0 ->
    0;
lbinom(N, K) when N == K ->
    0;
lbinom(N, K) ->
    kf_lgamma(N + 1) - kf_lgamma(K + 1) - kf_lgamma(N - K + 1).

%% Hypergeometric distribution
%%
%%  n11  n12  | n1_
%%  n21  n22  | n2_
%% -----------+----
%%  n_1  n_2  | n
%%
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

%% Calculate fisher exact test
%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Helper functions used by fisher_exact
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Helper function to sum up probabilities for the left tail.
sumleft(Left, I, Aux, P, Q, Max) when P < 0.99999999 * Q, I =< Max ->
    {Aux1, P1} = hypergeo_acc(I, 0, 0, 0, Aux, P),
    sumleft(Left + P, I + 1, Aux1, P1, Q, Max);
sumleft(Left, I, Aux, P, _, _) ->
    {Left, P, I, Aux}.

%% Helper function to sum up probabilities for the right tail.
sumright(Right, J, Aux, P, Q) when P < Q, J >= 0 ->
    {Aux1, P1} = hypergeo_acc(J, 0, 0, 0, Aux, P),
    sumright(Right + P, J - 1, Aux1, P1, Q);
sumright(Right, J, Aux, P, _) ->
    {Right, P, J, Aux}.

%% Helper function to adjust the left tail.
adjustleft(P, Q, Left, I) when P < 1.00000001 * Q ->
    {Left + P, I};
adjustleft(_, _, Left, I) ->
    {Left, I - 1}.

%% Helper function to adjust the right tail.
adjustright(P, Q, Right, J) when P < 1.00000001 * Q ->
    {Right + P, J};
adjustright(_, _, Right, J) ->
    {Right, J + 1}.

%% Calculate the twosided statistic.
gettwo(Left, Right) when Left + Right >= 1.0 ->
    1.0;
gettwo(Left, Right) ->
    Left + Right.

%% Adjust left and right
adjust_left_right(Left, _, I, J, N11, Q) when abs(I - N11) < abs(J - N11) ->
    Right1 = 1.0 - Left + Q,
    {Left, Right1};
adjust_left_right(_Left, Right, _, _, _, Q) ->
    {1.0 - Right + Q, Right}.


