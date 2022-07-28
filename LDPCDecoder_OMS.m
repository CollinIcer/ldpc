function [ v ] = LDPCDecoder_OMS( H, LLR_y, beta, iterMax,mb )
    %LDPCDecoder_SP LDPC使用积和算法解码
    %   H为校验矩阵，用以进行判决；LLR_y为接收到的信号初始置信度；iterMax为最大迭代次数；返回v为解码后的信息序列估计�?
    %初始�?u �?v 矩阵
    U0i = LLR_y;
    Uji = zeros(size(H));
    Vij = zeros(size(H'));
    [VerificationNodes, VariableNodes] = size(H);
    x = zeros(size(LLR_y));
    
    for iter = 1:1:iterMax
       %disp(['the ' num2str(iter) '-th iteration of OMS'])
        % 求解Vij矩阵
        for i = 1:1:VariableNodes
            idx = find(H(:, i) == 1);
            for k = 1:1:length(idx)
                 Vij(i, idx(k)) =  U0i(i) + sum(Uji(idx, i)) - Uji(idx(k), i);
            end
        end
        
        % 求解Uji矩阵
        for j = 1:1:VerificationNodes
            idx = find(H(j, :) == 1);
            for k = 1:1:length(idx)
                %Uji(j, idx(k)) = 2*atanh(prod(tanh(Vij(idx, j)/2))/tanh(Vij(idx(k), j)/2));
                % prodOfSign = prod( sign(Vij(idx, j)) ) / sign(Vij(idx(k), j));
                 % 排除掉Vij(idx(k), j)，求符号
                prodVal = 1.0;
                for t = 1:1:length(idx)
                    if t ~= k
                        prodVal = prodVal * sign( Vij(idx(t), j) );
                    end
                end
                % 排除掉Vij(idx(k), j)，求min
                if k == 1
                    minOfVal = min(abs(Vij(idx(2:end), j)));
                elseif k == length(idx)
                    minOfVal = min(abs(Vij(idx(1:k-1), j)));
                else
                    minOfVal = min(min( abs(Vij(idx(1:k-1), j)) ), min( abs(Vij(idx(k+1:end), j)) ) );
                end
                % Uji(j, idx(k))赋�?
                Uji(j, idx(k)) = prodVal * max(minOfVal - beta, 0);
            end
        end
        
        %判决
        for i = 1:1:length(x)
            idx = find(H(:, i) == 1);
            addVal = sum(Uji(idx, i)) + U0i(i);
            if(addVal < 0)
                x(i) = 1;
            else
                x(i) = 0;
            end
        end
        
        %如果校验关系满足 break;
        %否则继续迭代
        if mod(H*(x'), 2) == 0
            break;
        end
    end
    
    v = x(mb+1:end);
    
    end
    
    