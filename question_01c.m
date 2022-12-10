%% Question 1 of assignment 2
% Shway Wang
% Dec. 8, 2022

%% Preliminaries
% the total weight
total_weight = 1;

% lower bound
p_min = 10;

% upper bound
p_max = 20;

% epsilon, weight's upper bound
epsilon = 0.01;

% number of items
n = 1000;

% number of rounds
numRounds = 1;

% the verbose
verbose = false;

% the best CR
bestCR = 1 + log(p_max / p_min);

totalER = 0;
for roundI = 1:numRounds
    % initialize the input sequence
    input_seq = init_input_seq(p_min, p_max, epsilon, n);
    v_list = zeros(n, 1);
    w1_list = zeros(1, n);
    w2_list = zeros(1, n);
    iI = 1;
    for vI = 1:3:length(input_seq)
        v_list(iI) = -1 * input_seq(vI);
        w1_list(iI) = input_seq(vI + 1);
        w2_list(iI) = input_seq(vI + 2);
        iI = iI + 1;
    end
    
    % compute the answer by the threshold-based algorithm
    decision_list = threshold_based_alg(p_min, p_max, input_seq, total_weight);

    % get optimal solution by solving LP
    c = v_list;
    intcon = 1:length(v_list);
    
    A = [w1_list; w2_list];
    b = [total_weight; total_weight];
    
    options = optimoptions('intlinprog');
    options.IntegerTolerance = 0.000001;
    x = intlinprog(c, intcon, A, b, [], [], zeros(1, n), ones(1, n), [], options);
    
    % display the computed results
    opt_v = 0;
    opt_w = 0;
    alg_v = 0;
    alg_w1 = 0;
    alg_w2 = 0;
    for i = 1:length(v_list)
	    opt_v = opt_v + -1 * v_list(i) * x(i);
	    opt_w = opt_w + w1_list(i) * x(i);
	    alg_v = alg_v + -1 * v_list(i) * decision_list(i);
	    alg_w1 = alg_w1 + w1_list(i) * decision_list(i);
        alg_w2 = alg_w2 + w2_list(i) * decision_list(i);
    end
    
    if (verbose)
        fprintf("opt_w = %0.7f; total_weight = %0.7f\n", opt_w, total_weight); %#ok<UNRCH> 
        fprintf("alg_w1 = %0.7f; total_weight = %0.7f\n", alg_w1, total_weight);
        fprintf("alg_w2 = %0.7f; total_weight = %0.7f\n", alg_w2, total_weight);
    end
    
    assert(opt_w <= total_weight + 0.0001);
    assert(alg_w1 <= total_weight);
    assert(alg_w2 <= total_weight);
    
    er = opt_v/alg_v;
    if (verbose)
        fprintf('opt: %0.7f\n', opt_v); %#ok<UNRCH> 
        fprintf('alg: %0.7f\n', alg_v);

        % emperical ratio = opt / alg
        fprintf('emperical ratio: %0.7f\n', er);
    end

    totalER = totalER + er;
end
fprintf('average emperical ratio: %0.7f\n', totalER / numRounds);
fprintf('best CR in theory: %0.7f\n', bestCR);
assert(opt_v/alg_v <= bestCR);


%%%%%%%%%%%%%%%%%%%%%%%% Aux Funcs %%%%%%%%%%%%%%%%%%%%%%%%%%%

function decision_seq = threshold_based_alg(p_min, p_max, input_seq, total_weight)
%% Implement the threshold-based algorithm
% Shway Wang
% Dec. 8, 2022

arguments
    p_min (1,1) double
    p_max (1,1) double
    input_seq (1,:) double
    total_weight (1,1) double
end

y1 = 0;
y2 = 0;
decision_seq = [];
for tI = 1:3:length(input_seq)
    v = input_seq(tI);
    w1 = input_seq(tI + 1);
    w2 = input_seq(tI + 2);
    p = threshold(p_min, p_max, y1, y2, total_weight);
    x = 0;
    if (v/w1 >= p && v/w2 >= p)
        x = 1;
    end

    if (y1 + w1 * x > total_weight || y2 + w2 * x > total_weight)
        x = 0;
    end

    y1 = y1 + w1 * x;
    y2 = y2 + w2 * x;
    decision_seq = [decision_seq, x]; %#ok<AGROW> 
end
end

function th = threshold(p_min, p_max, y1, y2, total_weight)
%% Compute the threshold value for given bounds
% Shway Wang
% Dec. 8, 2022

arguments
    p_min (1,1) double
    p_max (1,1) double
    y1 (1,1) double
    y2 (1,1) double
    total_weight (1,1) double
end

beta = 1 / (1 + log(p_max/p_min));
if (y1 >= 0 && y1 < beta && y2 >= 0 && y2 < beta)
    th = p_min;
elseif ((y1 >= beta && y1 <= total_weight) || (y2 >= beta && y2 <= total_weight))
    th = p_min * exp(((y1 + y2) / 2) / beta - 1);
end
end

function seq = init_input_seq(p_min, p_max, epsilon, n)
%% Initializes the input sequence
% Shway Wang
% Dec. 8, 2022

arguments
    p_min (1,1) double
    p_max (1,1) double
    epsilon (1,1) double
    n (1,1) double
end

% initialize empty sequence
seq = zeros(1, n * 3);
for i = 1:n
    w = unifrnd(0, epsilon);
    while (w == 0)
        w = unifrnd(0, epsilon);
    end

    w2 = unifrnd(0, epsilon);
    while (w2 == 0)
        w2 = unifrnd(0, epsilon);
    end
    
    % initialize v randomly
    ratio = unifrnd(p_min, p_max);
    v = w * ratio;
    
    assert(v/w >= p_min && v/w <= p_max);
    
    % append v and w to seq
    seq((i - 1) * 3 + 1) = v;
    seq((i - 1) * 3 + 2) = w;
    seq((i - 1) * 3 + 3) = w2;
end
end


