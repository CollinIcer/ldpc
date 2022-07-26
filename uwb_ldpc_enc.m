function [ x ] = Encoder_uwb_ldpc(s)


    %-------------------------------------Hs, Hp----------
    %% Pre Definition
    n = 2304; % 编码序列总比特数 or H矩阵列数
    k = 1920; % 信息序列比特数 
    r = k/n; % 编码速率
    z = 96; % Hij_Matrix方阵大小
    mb= (2304-1920)/z;
    nb=n/z;
    kb=nb-mb;
    %% Coding Algorithm 2
    % Generate H, Hp, Hs
    H = zeros(n-k, n);
    Hij_Matrix = zeros(z, z);
    
    H_Block=[ [1  25 55 -1 47 4  -1 91 84 8  86 52 82 33 5  0  36 20 4  77 80 0   -1 -1]
              [-1 6  -1 36 40 47 12 79 47 -1 41 21 12 71 14 72 0  44 49 0  0  0   0  -1]
              [51 81 83 4  67 -1 21 -1 31 24 91 61 81 9  86 78 60 88 67 15 -1 -1  0   0]
              [68 -1 50 15 -1 36 13 10 11 20 53 90 29 92 57 30 84 92 11 66 80 -1 -1   0]];
    
    
    % Generate H
    for i = 1:1:mb
        for j = 1:1:nb
           idx = H_Block(i, j)+1;
           Hij_Matrix = zeros(z, z);
           if(idx == 0) %XKE FIX
               continue;
           end
           Hij_Matrix(1, idx) = 1;
           %if(idx == 56 && j <= kb)
           %    Hij_Matrix(1, idx) = 0;
           %end
           
           for k = 2:1:z
               idx = idx + 1;
               if(idx > z)
                   idx = 1;
               end
               Hij_Matrix(k, idx) = 1;
           end
           H((i-1)*z+1:i*z, (j-1)*z+1:j*z) = Hij_Matrix;
        end
    end
    
    % Generate Hp, Hs
    % Hp = H(1:mb*z, 1:mb*z);
    % Hs = H(1:mb*z, mb*z+1:nb*z);
    Hs = H(1:mb*z, 1:kb*z);
    Hp = H(1:mb*z,kb*z+1:nb*z );
    
    
    %---------------------------------------------------
    
    
    p = zeros(mb*z,1); %校验比特序列p
    
    if 1
        sw=zeros(z,1);
        sw_cmp = zeros(z,mb);
        sw_tmp = zeros(z,mb);
        for i=1:1:kb
            for j=1:1:mb
             wj=Hs((j-1)*z+1:j*z, (i-1)*z+1:i*z);
             info_block = s((i-1)*z+1:i*z);
             swj = wj * info_block';
             sw = mod((swj + sw),2);
             sw_tmp(1:z, j) = mod(swj ,2);
            end
            sw_cmp= mod((sw_tmp +  sw_cmp),2);
        end
        
       
        
        %sum each  of sw
        p1 = zeros(z,1);
        for j=1:1:z  
            for i = 1:1:mb
                p1(j,1)= mod((sw_cmp(j,i)+p1(j,1)),2);
            end
        end
        p(1:z,1)=p1;
        
        if p1==sw
            display("P1 == sw_cmp");
        end
        
        for i=2:1:mb  
            srow=zeros(z,1);
            for j=1:1:kb
                sel_Hs = Hs((i-2)*z+1:(i-1)*z,(j-1)*z+1:j*z);
                sel_s = s((j-1)*z+1:j*z);
                srow = mod((sel_Hs*sel_s')+srow,2);
            end
            if(srow == sw_cmp(1:z,i-1))
                display("srow already cal");
            end
            
            if(i==2)
                p1_mul = Hp((i-2)*z+1:(i-1)*z,1:z) * p(1:z,1);
                p((i-1)*z+1:i*z,1) = mod((srow+p1_mul),2);
            else if(i==3)%r+1, r=2
                p((i-1)*z+1:i*z,1) = mod((srow+p((i-2)*z+1:(i-1)*z,1) + p(1:z,1)),2);
              
            else
                 p((i-1)*z+1:i*z,1) = mod((srow+p((i-2)*z+1:(i-1)*z,1)),2);
            end
           end
        end
         
        
        
    end
    
    
    
    
    x = [p' s];
    %x序列的正确性检验
    H= [Hp Hs];
    y=mod(H*x',2);
    t=sum(y);
    if t ~= 0
        display("ERROR in ENC");
    end
        
        
    
    end